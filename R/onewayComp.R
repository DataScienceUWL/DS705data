onewayComp <- function(formula,data,alpha=.05,var.equal=TRUE,con=NA,nboot=0,adjust='one.step'){
  #
  # formula = x~g  is the one-way formula, x is response, g is factor variable (for example)
  # data = data frame with x and g (or whatever you've name the response and predictor)
  # alpha = family-wise error rate (conf = 1 - alpha)
  # var.equal = TRUE for pooled variance (one.step -> Tukey-Kramer) or FALSE for unpooled (one.step -> Games-Howell)
  # con = contrast matrix, if unspecified or NA, then all pairs are assumed, otherwise 
  #       con is (number of contrasts) by (number of groups) and each row is a contrast
  #       rownames from con will be displayed if available
  # nboot = number of bootstrap resamples (null restricted bootstrap)
  # adjust = 'one.step' (default), 'bonferroni', 'holm', 'none' (or anything else that can be passed to p.adjust() )
    
  if (missing(data)) {
    xdf <- model.frame(formula)
  } else {
    xdf <- model.frame(formula, data)
  }
  cl <- match.call()
  # note model.frame will remove NA by default
  
  # make sure group labels are a factor
  xdf[,2] <- factor(xdf[,2])
  levs <- levels(xdf[,2])
  J <- length(levs)
  
  mean.obs <- tapply(xdf[,1],xdf[,2],mean)
  n.obs <- tapply(xdf[,1],xdf[,2],length)
  v.obs <- tapply(xdf[,1],xdf[,2],var)
  se2.obs <- v.obs/n.obs
  N <- length(xdf[,1])
  pooledVar.obs <- sum((v.obs*(n.obs-1)))/(N-J) # pooled variance estimate
  
  # center the samples to satisfy null hypothesis (we are bootstrapping the residuals)
  xlist.ctr <- split(xdf[,1],xdf[,2])
  for (i in 1:J){
    xlist.ctr[[i]] <- xlist.ctr[[i]]-mean.obs[i]
  }
  
  if (anyNA(con)){
    customCon = FALSE
    # build pairwise contrast matrix
    num.compare <- J*(J-1)/2
    con <- matrix(0,num.compare,J)
    rowLabel <- character(num.compare)
    ind <- 1
    for (i in 1:(J-1)){
      for (j in (i+1):J){
        con[ind,i] <- -1
        con[ind,j] <- 1
        rowLabel[ind] <- paste0(levs[j],'-',levs[i])
        ind <- ind + 1
      }
    }
  } else {
    customCon = TRUE
    num.compare <- as.numeric(dim(con)[1])
    rowLabel <- rownames(con)
    if (is.null(rowLabel)){
      rowLabel <- rep("",num.compare)
    } 
    for (j in 1:num.compare){
      if (rowLabel[j]==""){
        rowLabel[j] <- paste('contrast',j)
      }
    }
  }
  
  if (customCon & adjust=='one.step'){
    print('****** Warning ******')
    print('you have given a custom contrasts matrix')
    print('and asked for Tukey or Games-Howell single step p-values and CIs')
    print('but Tukey and GH are only for all pairwise comparisons ...')
    print('so the results may be nonsense.')
    print('If you want to use Tukey or GH for all pairs, specify con=NA')
    print('*********************')
  }
  
  if (nboot > 0){
    # resampling and collecting stats
    # sample residuals with replacment and within each group
    mean.samp <- matrix(0,nboot,J)
    v.samp <- matrix(0,nboot,J)
    
    for (i in 1:J){
      x <- matrix(sample(xlist.ctr[[i]],n.obs[i]*nboot,replace=TRUE),nboot,n.obs[i])
      mean.samp[,i] <- apply(x,1,mean)
      v.samp[,i] <- apply(x,1,var)
    }
    x <- NULL
    
    # now build pairwise bootstrapped t-test statistics
    tstat <- matrix(0,nboot,num.compare)
    if (var.equal==TRUE){ # computed pooled variance estimate
      pooledVar.samp <- (v.samp %*% (n.obs-1))/(N-J)  # this is MSE
    }
    
    for (i in 1:num.compare){
      w <- con[i,] # get weights for contrast i
      if (var.equal==TRUE){
        tstat[,i] <- (mean.samp %*% w)/sqrt( pooledVar.samp*(sum((w^2)/n.obs)))
      } else {
        tstat[,i] <- (mean.samp %*% w)/sqrt( v.samp %*% ((w^2)/n.obs))
      }
    }
    
    # now construct Studentized range distribution
    qstat <- apply(abs(tstat),1,max)
    
    # get critical value
    qcrit <- quantile(qstat,p=1-alpha)
    
    # now build final intervals
    ci <- matrix(0,num.compare,2)
    psihat <- numeric(num.compare)
    padjust <- numeric(num.compare)
    p <- numeric(num.compare)
    t.obs <- numeric(num.compare)
    for (j in 1:num.compare){
      w <- con[j,]
      mean.diff <- mean.obs %*% w
      psihat[j] <- mean.diff
      if (var.equal==TRUE){
        se <- sqrt( pooledVar.obs*(sum((w^2)/n.obs)))
      } else {
        se <- sqrt( sum((v.obs*(w^2))/n.obs) )
      }
      
      if (adjust=='one.step'){
        ci[j,1] <- mean.diff - qcrit*se
        ci[j,2] <- mean.diff + qcrit*se
      } else if (adjust=='bonferroni') {
        alphaP2 <- alpha/num.compare/2
        tcrit <- unname(quantile(tstat[,j],p=c(alphaP2,1-alphaP2)))
        ci[j,1] <- mean.diff - tcrit[2]*se
        ci[j,2] <- mean.diff - tcrit[1]*se
      } else if (adjust=='none') {
        alpha2 <- alpha/2
        tcrit <- unname(quantile(tstat[,j],p=c(alpha2,1-alpha2)))
        ci[j,1] <- mean.diff - tcrit[2]*se
        ci[j,2] <- mean.diff - tcrit[1]*se
      } else {
        ci[j,1] <- NA
        ci[j,2] <- NA
      }
      t.obs[j] <- as.numeric(mean.diff/se)
      p[j] <- sum( abs(tstat[,j]) > abs(t.obs[j]) )/nboot
      padjust[j] <- sum( qstat>abs(t.obs[j]) )/nboot
    }
  } else { # not bootstrapping
    # now build final intervals
    reduceDF = FALSE # FALSE gives pooled df for errors when var.equal=TRUE
    ci <- matrix(0,num.compare,2)
    psihat <- numeric(num.compare)
    padjust <- numeric(num.compare)
    p <- numeric(num.compare)
    t.obs <- numeric(num.compare)
    for (j in 1:num.compare){
      w <- con[j,]
      mean.diff <- mean.obs %*% w
      psihat[j] <- mean.diff
      if (var.equal==TRUE){
        se <- sqrt( pooledVar.obs*(sum((w^2)/n.obs)))
        ind <- w!=0
        if ( reduceDF & (sum(ind)==2) ){
          t.df <- sum(n.obs[ind])-2
        } else {
          t.df <- N-J
        }
      } else {
        se <- sqrt( sum((v.obs*(w^2))/n.obs) )
        # compute welch correct df for Games-Howell or t intervals
        t.df <- sum(se2.obs*abs(w))^2/( sum(se2.obs^2*abs(w)/(n.obs-1)))
      }
      t.obs[j] <- as.numeric(mean.diff/se)
      p[j] <- 2*pt(-abs(t.obs[j]),t.df)
      if (adjust=='one.step'){
        # one.step only makes sense for pairwise comparisons
        qcrit <- qtukey(1-alpha,J,t.df)/sqrt(2)
        ci[j,1] <- mean.diff - qcrit*se
        ci[j,2] <- mean.diff + qcrit*se
        padjust[j] <- ptukey(sqrt(2)*abs(t.obs[j]),lower.tail=FALSE,df=N-J,nmeans=J)
      } else if (adjust=='bonferroni') {
        alphaP2 <- alpha/num.compare/2
        tcrit <- qt(1-alphaP2,t.df)
        ci[j,1] <- mean.diff - tcrit*se
        ci[j,2] <- mean.diff + tcrit*se
      } else if (adjust=='none') {
        alpha2 <- alpha/2
        tcrit <- qt(1-alpha2,t.df)
        ci[j,1] <- mean.diff - tcrit*se
        ci[j,2] <- mean.diff + tcrit*se
      } else {
        ci[j,1] <- NA
        ci[j,2] <- NA
      }
    }
  }

  if (adjust=='holm'){
    padjust <- p.adjust(p,'holm')
    reject <- (padjust<alpha)
    comp.matrix <- cbind(psihat,ci,t.obs,p,padjust,reject)
    dimnames(comp.matrix) <- list(rowLabel,c('diff','lwr','upr','t','p','p adj','rej H_0'))
  } else if (adjust=='one.step'){
    reject <- padjust < alpha
    comp.matrix <- cbind(psihat,ci,t.obs,p,padjust,reject)
    dimnames(comp.matrix) <- list(rowLabel,c('diff','lwr','upr','t','p','p adj','rej H_0'))
  } else if (adjust=='bonferroni'){
   padjust <- p.adjust(p,adjust)
   reject <- padjust < alpha
   comp.matrix <- cbind(psihat,ci,t.obs,p,padjust,reject)
   dimnames(comp.matrix) <- list(rowLabel,c('diff','lwr','upr','t','p','p adj','rej H_0'))
  } else if (adjust=='none'){
    padjust <- p
    reject <- p < alpha
    comp.matrix <- cbind(psihat,ci,t.obs,p,reject)
    dimnames(comp.matrix) <- list(rowLabel,c("diff","lwr","upr","t","p",'rej H_0'))
  } else {
    padjust <- p.adjust(p,adjust)
    comp.matrix <- cbind(psihat,ci,t.obs,p,padjust)
    dimnames(comp.matrix) <- list(rowLabel,c("diff","lwr","upr","t","p","p adj"))
  }
  
  if (!any(is.na(con))){
    # now build a pairwise matrix like output from pairwise.t.test
    ind <- 1
    p.value <- matrix(NA,J-1,J-1)
    for (i in 1:(J-1)){
      for (j in i:(J-1)){
        p.value[j,i] <- padjust[ind]
        ind <- ind + 1
      }
    }
    dimnames(p.value) <- list(levs[2:J],levs[1:(J-1)])
    p.adjust.method <- adjust
    if (var.equal){
      method <- 't tests with pooled SD'
    } else {
      method <- 't tests with unpooled SD'
    }
    nms <- names(xdf)
    data.name <- paste(nms[1],"and",nms[2])
    pairwt <- list(method=method,data.name=data.name,
                   p.value=p.value,p.adjust.method=p.adjust.method)
    class(pairwt) <- 'pairwise.htest'
    
    result <- list(call=cl,comp=comp.matrix,pairw=pairwt)
  } else {
    result <- list(call=cl,comp=comp.matrix)
  }
  result
}
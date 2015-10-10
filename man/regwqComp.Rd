\name{regwqComp}
\alias{regwqComp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
reqwqComp
}
\description{
Implements the multistep REGWQ procedure described in the DS705 presentations.
}
\usage{
regwqComp(formula, data, alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
%%     ~~Describe \code{formula} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (formula, data, alpha = 0.05) 
{
    if (missing(data)) {
        dat <- model.frame(formula)
    }
    else {
        dat <- model.frame(formula, data)
    }
    if (ncol(dat) != 2) {
        stop("Specify one response and only one class variable in the formula")
    }
    if (is.numeric(dat[, 1]) == FALSE) {
        stop("Response variable must be numeric")
    }
    response <- dat[, 1]
    group <- as.factor(dat[, 2])
    fl <- levels(group)
    a <- nlevels(group)
    N <- length(response)
    samples <- split(response, group)
    n <- sapply(samples, "length")
    mm <- sapply(samples, "mean")
    vv <- sapply(samples, "var")
    MSE <- sum((n - 1) * vv)/(N - a)
    df <- N - a
    nc <- a * (a - 1)/2
    order.h1 <- data.frame(Sample = fl, SampleNum = 1:a, Size = n, 
        Means = mm, Variance = vv)
    ordered <- order.h1[order(order.h1$Means, decreasing = FALSE), 
        ]
    rownames(ordered) <- 1:a
    i <- 1:(a - 1)
    h1 <- list()
    for (s in 1:(a - 1)) {
        h1[[s]] <- i[1:s]
    }
    vi <- unlist(h1)
    j <- a:2
    h2 <- list()
    for (s in 1:(a - 1)) {
        h2[[s]] <- j[s:1]
    }
    vj <- unlist(h2)
    h3 <- list()
    h4 <- list()
    for (s in 1:(a - 1)) {
        h3[[s]] <- rep(j[s], s)
        h4[[s]] <- rep(i[s], s)
    }
    Nmean <- unlist(h3)
    Step <- unlist(h4)
    mean.difference <- sapply(1:nc, function(arg) {
        i <- vi[arg]
        j <- vj[arg]
        (ordered$Means[j] - ordered$Means[i])
    })
    T <- sapply(1:nc, function(arg) {
        i <- vi[arg]
        j <- vj[arg]
        (ordered$Means[j] - ordered$Means[i])/sqrt(MSE * (1/ordered$Size[i] + 
            1/ordered$Size[j]))
    })
    pvalues <- ptukey(T * sqrt(2), Nmean, df, lower.tail = FALSE)
    alpha.level <- 1 - (1 - alpha)^(Nmean/a)
    level1 <- (Nmean == a)
    level2 <- (Nmean == a - 1)
    level3 <- level1 + level2
    alpha.level[level3 == 1] <- alpha
    quantiles <- qtukey(1 - alpha.level, Nmean, df)
    for (h in 1:(nc - 1)) {
        if (quantiles[h + 1] >= quantiles[h]) {
            quantiles[h + 1] <- quantiles[h]
        }
    }
    Rejected1 <- (pvalues < alpha.level)
    names.ordered <- sapply(1:nc, function(arg) {
        i <- vi[arg]
        j <- vj[arg]
        paste(ordered$Sample[j], "-", ordered$Sample[i], sep = "")
    })
    for (s in 1:nc) {
        if (Rejected1[s] == FALSE) {
            Under1 <- (vj[s] >= vj)
            Under2 <- (vi[s] <= vi)
            Under3 <- Under1 * Under2
            Under4 <- which(Under3 == 1)
            Rejected1[Under4] <- FALSE
        }
    }
    rowOrderMat <- matrix(0, a, a)
    rowOrderMat[lower.tri(rowOrderMat)] <- 1:nc
    rownames(rowOrderMat) <- fl
    colnames(rowOrderMat) <- fl
    rowOrderVec <- numeric(nc)
    signVec <- rep(1, nc)
    for (s in 1:nc) {
        i <- vi[s]
        j <- vj[s]
        si <- ordered$SampleNum[i]
        sj <- ordered$SampleNum[j]
        if (si < sj) {
            rowOrderVec[s] <- rowOrderMat[sj, si]
        }
        else {
            rowOrderVec[s] <- rowOrderMat[si, sj]
            names.ordered[s] <- paste(ordered$Sample[i], "-", 
                ordered$Sample[j], sep = "")
            mean.difference[s] <- -mean.difference[s]
            T[s] <- -T[s]
        }
    }
    ind <- order(rowOrderVec)
    pvalues <- pvalues[ind]
    mean.difference <- mean.difference[ind]
    T <- T[ind]
    names.ordered <- names.ordered[ind]
    Rejected1 <- Rejected1[ind]
    alpha.level <- alpha.level[ind]
    pv <- 2 * pt(-abs(T), df = N - a)
    comp.matrix <- cbind(mean.difference, T, pv, pvalues, alpha.level, 
        as.numeric(Rejected1))
    dimnames(comp.matrix) <- list(names.ordered, c("diff", "t", 
        "p", "p adj", "alpha adj", "rej H_0"))
    return(comp.matrix)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

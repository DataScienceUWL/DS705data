\name{onewayComp}
\alias{onewayComp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simultaneous Confidence Intervals and Hypothesis Tests for One-way Designs
}
\description{
Generates Studentized confidence intervals and/or hypothesis tests for one-way designs.  P-values can be adjusted ...
}
\usage{
onewayComp(formula,data,alpha=.05,var.equal=TRUE,con=NA,nboot=0,adjust='one.step')
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{e.g. x~group, this would the same formula you would supply to lm for a one-way ANOVA}
  \item{data}{data frame containing x and group}
  \item{alpha}{significance level of the test, 1-alpha is the overall confidence level}
  \item{var.equal}{Defaults to TRUE, if FALSE, then Welch type corrections are used for the Studentized test statistics.  FALSE along with adjust='one.step' gives the Games-Howell procedure.}
  \item{con}{optional contrast matrix.  If con=NA, then all pairwise comparisons are made.}
  \item{nboot}{if nboot=0, then t-distributions are used for P-values and critical values (including Tukey max t), but if nboot>0 then resampling is used to approximate the critical values and P-values.  For the bootstrap resampling separate samples are used and each sample is recentered so that the null hypothesis is true.}
  \item{adjust}{possible values are 'one.step' for Tukey Kramer or Games Howell procedures that give both hypothesis tests and confidence intervals (cannot be used with custom contrast matrix). 'bonferroni', 'none', 'holm', or any other adjustment that can be passed to p.adjust()}
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{call}{Shows the details of the function call.}
\item{comp}{A matrix, with names columns.  'diff' is the sample mean difference or sample mean contrast value. 'lwr' and 'upr' are the confidence interval limits (if available). 't' is the value of the Studentized test statistic. 'p' is the unadjusted p-value as would be computed by a single t-test. 'p adj' is the adjusted p-value as determined by the selected adjust method.  'Rej H_0' is 1 if p adj is less than alpha}
\item{pairw}{Matrix showing pairwise comparisons and P-values.}
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
Jeff Baggett, University of Wisconsin - La Crosse.
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
function (x) 
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

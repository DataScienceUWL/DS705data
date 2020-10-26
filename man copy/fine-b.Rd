\name{fine-b}
\alias{fine-b}
\docType{data}
\title{
Weight gain
%%   ~~ data name/kind ... ~~
}
\description{Weight gain data for rats exposed to aconiazide, 
 a proposed treatment for tuberculosis.
   Taken from Fine and Bosch, JASA, to appear 2000

   For details of the study which generated these data, see
    Dooley et al., Journal of the American College of Toxicology, 1993
   For additional analyses of these data, see
    Kodell and West, Risk Analysis, 1993    
    Bosch et al., ASA Proceedings of the Biometrics Section, 1995
    Bosch et al., Risk Analysis, 1996
   
   Thanks to Frederick Beland for making these data available.


 W is weight gain over the 14-day duration  
  of the experiment, on 50 female Fischer 344 rats 
 DOSE is the dose of aconiazide (mg/kg/day)
}
\usage{data("fine-b")}
\format{
  A data frame with 50 observations on the following 2 variables.
  \describe{
    \item{\code{W}}{Weight gain.}
    \item{\code{Dose}}{Dose of aconiazide (mg/kg/day.}
      }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(fine-b)
## maybe str(fine-b) ; plot(fine-b) ...
}
\keyword{datasets}

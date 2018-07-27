\name{rowtime}
\alias{rowtime}
\docType{data}
\title{
Row Times 
%%   ~~ data name/kind ... ~~
}
\description{Ninety members (aged = 18.1 – 23.4 years) of three Division I women’s intercollegiate rowing teams (National Collegiate Athletic Association) within the Big Ten Conference volunteered to participate in a study to predict race time for female collegiate rowers from physical characteristics.  Some of the measured characteristics have been omitted from this set.
}
\usage{data("rowtime")}
\format{
  A data frame with 20 observations on the following 2 variables.
  \describe{
    \item{\code{racetime}}{2K rowing ergometer time, sec.}
    \item{\code{tall}}{Height, cm.}
    \item{\code{weight}}{Weight, kg.}
    \item{\code{arm_span}}{Arm span, cm.}
    \item{\code{flex_arm}}{Arm muscle circumference, cm.}
    \item{\code{thigh_ci}}{Thigh circumference, cm}
    \item{\code{calf_cir}}{Calf circumference, cm}
    \item{\code{tricep}}{Triceps skinfold, mm}
    \item{\code{biceps}}{Biceps skinfold, mm}
    \item{\code{thigh}}{Thigh skinfold, mm}
    \item{\code{est_ffm}}{Estimated fat-free mass, kg}
    \item{\code{est_fm}}{Estimated fat mass, kg}
    \item{\code{best_snr}}{Best sit-n-reach, cm}
    \item{\code{best_vj}}{Best vertical jump, cm}
    \item{\code{leg_powe}}{Le power, kpm/min}
    \item{\code{endo}}{Endomorphy}
    \item{\code{meso}}{Mesomorphy}
    \item{\code{ecto}}{Ectomorphy}
    \item{\code{exp_var}}{Collegiate experience, 1=Yes, 0=No}
    \item{\code{pre_exp}}{Pre-collegiate experience, 1=Yes, 0=No}
      }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{Based on "Use of Modeling to Predict Performance from Physical Characteristics in Collegiate Female Rowers" by Battista, et.al.
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(rowtime)
## maybe str(rowtime) ; plot(rowtime) ...
}
\keyword{datasets}

\name{HealthExam}
\alias{HealthExam}
\docType{data}
\title{
HealthExam data set.
%%   ~~ data name/kind ... ~~
}
\description{
A dataset containing weights, blood pressures and other attributes of 40 adult men and 40 adult women.  Sample is loosely based on data in the National Health Exam III.
}
\usage{data("HealthExam")}
\format{
  A data frame with 80 observations on the following 12 variables.
  \describe{
    \item{\code{Sex}}{male or female}
    \item{\code{AgeGroup}}{18-35, 36-65, 65+}
    \item{\code{Height}}{height in inches}
    \item{\code{Weight}}{weight in pounds}
    \item{\code{Waist}}{circumference in cm}
    \item{\code{Pulse}}{pulse rate in beats per minute}
    \item{\code{SysBP}}{systolic blood pressure in mmHg}
    \item{\code{DiasBP}}{diastolic blood pressure in mmHg}
    \item{\code{BodyMass}}{body mass index in kg/m^2}
    \item{\code{Cholesterol}}{in mg/dL}
    \item{\code{Leg}}{upper leg length in cm}
    \item{\code{Region}}{Geographic region in the US.  1: Midwest, 2: Northeast, 3: South, 4: West}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{based on 
http://www.cdc.gov/nchs/nhanes/nh3data.htm
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(HealthData)
}
\keyword{datasets}

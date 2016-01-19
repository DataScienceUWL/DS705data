\name{Cryotherapy}
\alias{Cryotherapy}
\docType{data}
\title{Cryotherapy data set.}
\description{
Data was collected for a sample of college students at a university in the Midwest.  One variable measured was agility both before and at several points in time afterward.  This data set contains the agility measure before cryotherapy and at 17 minutes after cryotherapy in a cold whirlpool. Agility is measured by the time to perform the following exercise: "At the tester’s signal, the subject sprinted forward 
9.14 M (10 yards) and touched the tip of the cone (B) 
with their right hand. Then they performed a lateral shuf-
fle to the left 4.57 m (five yards) and touched the tip of 
the cone (C) with the left hand. Subjects then continued to 
shuffle 9.14 m (10 yards) to the right and touched the tip 
of the cone (D) with their right hand. They then shuffled 
4.57 m (five yards) to the left and touched point B with 
their left hand. Finally, subjects back peddled 9.14 m (10 
yards), passing through the finish at point A. Time was 
measured using a hand-held stopwatch." (quote taken from the article referenced below)  
}
\usage{data("Cryotherapy")
}
\format{
  A data frame with 18 observations on the following 2 variables.
  \describe{
    \item{\code{Precryotherapy}}{Agility measure prior to cryotherapy}
    \item{\code{After17Min}}{Agility measure 17 minutes after cryotherapy}
}}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\examples{
data(Cryotherapy)
}
\keyword{datasets}

\source{An author of "The effects of cold whirlpool on powe
r, speed, agility, and range of motion." Journal of Sports Science and Medicine (2008) 
7, 387-394}

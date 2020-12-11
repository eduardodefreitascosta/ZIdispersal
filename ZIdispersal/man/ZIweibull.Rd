\name{Ziweibull}
\alias{ZIweibull}
\title{ Ziweibull
%%   ~~ data name/kind ... ~~
}
\description{
This function runs a ZIweibull regression considering zeros and censored data
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{ZIweibull(t=dist,d=delta,data=wild_boar,co=c("age","sex"))}
\format{
  A data frame with none or 2 covariables.
  \describe{
    \item{\code{t}}{Dispersal distances}
    \item{\code{d}}{Indicative variable for censored observations}
    \item{\code{data}}{A dataset}
    \item{\code{co}}{A sting with the name sof covariables}
    \item{\code{dist}}{distribution to be used}
  }
}
\details{The user only can use 2 or none quantitative or binary covariables in this function
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\author{
Eduardo de Freitas Costa, Silvana Schneider, Giulia Carlotto, Taina Cabalheiro, Mauro Oliveira
%%  ~~who you are~~
}
\examples{
data(wild boar)
## maybe str(wild boar) ; plot(wild boar) ...
}
\keyword{datasets}

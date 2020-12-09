\name{ZIgamma}
\alias{ZIgamma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ZIgamma
%%  ~~function to do ... ~~
}
\description{This is a package to run a zero-inflated-censored dispersal distances regrssion. It uses a logistic and a gamma regressions to accomodate zeros and cenrored dispersal distances.
  Regression
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ZIgamma(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ze}{Indicative variabel for zero dispersal
  \item{t}{Dispersal distance variable
  \item{d}{Indicative variable for censured dispersal
  \item{data}{Dataset to be analysed
  \item{co}{ It is a string vector containg the names of covariates to be included in the model
%%     ~~Describe \code{x} here~~
}
}
}
}
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
 \item{flexsurvreg}
 \item{Surv}
 \item{glm}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Codes for the dispersal of wild boar ----
##-- ==>  read the data wild_boar,

## Run
ZIgamma(ze=zero,t=dist,d=delta,data=wild_boar,co=c("age","sex"))

{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.

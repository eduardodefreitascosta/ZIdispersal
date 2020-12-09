\name{ZIgamma}
\alias{ZIgamma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Zero inflated censored regression
%%  ~~function to do ... ~~
}
\description{This is a package to run a zero-inflated-censored dispersal distances regrssion. It uses a logistic and a gamma regressions to accomodate zeros and cenrored dispersal distances.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ZIgamma(ze,t,d,data,co=c())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ze}{Indicative variabel for zero dispersal}
  \item{t}{Dispersal distance variable}
  \item{d}{Indicative variable for censured dispersal}
  \item{data}{Dataset to be analysed}
  \item{co}{ It is a string vector containg the names of covariates to be included in the model}
%%     ~~Describe \code{x} here~~
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
\author{Eduardo Costa, Silvana Schneider, Giulia Carloto, Tainá Cabalheiro, Mauro Oliveira
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \item{}

%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Codes for the dispersal of wild boar ----
##-- ==>  read the data wild_boar, available on: https://github.com/eduardodefreitascosta/Ziweibull/blob/master/ZIgamma/wild_boar.rds

## Run
exemple<-ZIgamma(ze=zero,t=dist,d=delta,data=wild_boar,co=c("age","sex"))

The output is a list: 1 the logistic part; 2 gamma part. The user can ask for the summary for the first part as summary(exemmple$Logistic)

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

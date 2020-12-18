\name{ZIgamma}
\alias{ZIgamma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ZIgamma
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ZIgamma(ze,t,d,data,co=c())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ze}{
    Indicative variable for zero dispersal
}
  \item{t}{
    Dispersal observed distances
}

  \item{d}{
    Indicative variable for censored observation
}

  \item{data}{
    Dataset with the data
}

  \item{co}{
    A string with the names of the covariables used
}


}
\details{
Name of the variables in argument "co" must be equal to the dataset
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{Logistic}{Output of the logistic model}
  \item{Gamma}{Output of the gamma model}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Eduardo de Freitas Costa, Silvana Schneider, Giulia Carlotto, Taina Cabalheiro, Mauro Oliveira
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
flexruvreg, glm
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
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.

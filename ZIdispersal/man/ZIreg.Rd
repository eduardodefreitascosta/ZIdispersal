\name{ZIreg}
\alias{ZIreg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ZIreg
%%  ~~function to do ... ~~
}
\description{
This function runs a ZIgamma or a ZIweibull regression model to datasets with continous outcome, zero inflated and censored data.
}
\usage{
ZIreg(formula,zero,censor,data,dist)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under ‘Details’.
}
\item{zero}{
A vector containing an indicative variable used to indicate zeros.
}
\item{censor}{
A vector containing an indicative variable used to indicate zeros.
}
\item{data}{
an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which ZIreg is called.
}
\item{dist}{
a description of the distribution to be used in the model. For ZIreg this is a character string. only possible options are "weibull" and "gamma"
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
ZIreg returns different objects dependent on the distribution used.For "gamma" it is the output from a logistic glm object and a list of class "flexsurvreg". For more information see ?glm and ?flexsurvreg.
For dist="weibull", it is a list:

\item{call}{A copy of the function call}
\item{LogLik }{Log-likelihood}
\item{N }{The non-transformed maximum likelihood estimates}
\item{Zeros }{The number of zero observations}
\item{Censored }{The number of censored observations}
\item{df }{Degrees of freedom}
\item{AIC }{Akaike's information criterion (-2*log likelihood + 2*number of estimated parameters)}
\item{AICc }{Corrected Akaike's information criterion }
\item{BIC }{Bayesian information criterion }

%% ...
}
\references{
Jackson, C. (2016). flexsurv: A Platform for Parametric Survival Modeling in R. Journal of Statistical Software, 70(8), 1-33. doi:10.18637/jss.v070.i08

Louzada F, Moreira FF, de Oliveira MR (2018) A zero-inflated non defaultrate regression model for credit scoring data. Communications in Statistics-Theory and Methods 47(12):3002–3021

Nobre AA, Carvalho MS, Griep RH, Fonseca MdJMd, Melo ECP, Santos IdS,Chor D (2017) Multinomial model and zero-inflated gamma model to studytime spent on leisure time physical activity: an example of elsa-brasil. Revista de saude publica 51:76
}
\author{
Eduardo de Freitas Costa, Silvana Schneider, Giulia Bagatini Carlotto, Taina Cabalheiro, Mauro Oliveira
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Exemple of wild boar dispersal ----
url for download: https://github.com/eduardodefreitascosta/ZIdispersal/blob/master/ZIdispersal/wild_boar.rds

Zireg(dist~sex+age,zero=zero,censor=delta,data=wild_boar,dist="weibull")
or
Zireg(dist~sex+age,zero=zero,censor=delta,data=wild_boar,dist="gamma")

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

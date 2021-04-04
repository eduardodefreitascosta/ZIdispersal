\name{ZIreg}
\alias{ZIreg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ZIreg
%%  ~~function to do ... ~~
}
\description{
This function is used to fit ZIgamma or a ZIweibull regression model to datasets with continous outcome, zero-inflated and censored data.
}
\usage{
ZIreg(formula,zero,censor,data,dist)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under "Details".
}
\item{zero}{
a vector containing an indicative variable used to indicate zeros.
}
\item{censor}{
a vector containing an indicative variable used to indicate censoring.
}
\item{data}{
a data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which ZIreg is called.
}
\item{dist}{
a description of the distribution to be used in the model. For ZIreg this is a character string. only possible options are "weibull" and "gamma"
}
}
\details{
When dist="gamma", the model is a hurdle model using a logistic regression for zeros and gamma for distances larger than zero (Jackson, 2016; Nobre et al., 2017). The gamma model, by default, covariates are placed on both parameters of the distribution ,in this case the "shape" and "rate". When dist="weibull", the mode is a zero-inflated-censored weibull model (Louzada et al., 2018). A typical predictor has the form response ~ terms where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response. Also indicators for zero and for censoring must be provided in specific arguments. A terms specification of the form first + second indicates all the terms in first together with all the terms in second with any duplicates removed. A specification of the form first:second indicates the set of terms obtained by taking the interactions of all terms in first with all terms in second. The specification first*second indicates the cross of first and second. This is the same as first + second + first:second. The terms in the formula will be re-ordered so that main effects come first, followed by the interactions, all second-order, all third-order and so on: to avoid this pass a terms object as the formula.

}

\value{
ZIreg returns different objects dependent on the distribution used. For "gamma" it is the output from a logistic glm object and a list of class "flexsurvreg" (Jackson, 2016). For more information see ?glm and ?flexsurvreg.
The function "summary" can be applied to each object of the output: "$Logistic or $Gamma". 

For dist="weibull", the output is a list:

\item{call}{A copy of the function call}
\item{LogLik }{Log-likelihood}
\item{N }{The number of observations}
\item{Zeros }{The number of zero observations}
\item{Censored }{The number of censored observations}
\item{df }{Degrees of freedom}
\item{AIC }{Akaike's information criterion (-2*log likelihood + 2*k), where k is the number of predictors}
\item{AICc }{Corrected Akaike's information criterion (AIC+(2k(k+1))/(n-k)),where k is the number of predictors and n is the number of observations }
\item{BIC }{Bayesian information criterion (k*ln(n)-2*log likelihood), where k is the number of predictors and n is the number of observations}

%% ...
}
\references{
Jackson, C. (2016). flexsurv: A Platform for Parametric Survival Modeling in R. Journal of Statistical Software, 70(8), 1-33. doi:10.18637/jss.v070.i08

Louzada F, Moreira FF, de Oliveira MR (2018) A zero-inflated non defaultrate regression model for credit scoring data. Communications in Statistics-Theory and Methods 47(12):3002-3021

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
url for the dataset download: 
https://github.com/eduardodefreitascosta/ZIdispersal/blob/master/ZIdispersal/wild_boar.rds

Functions to run the model:
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

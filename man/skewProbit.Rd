\name{skewProbit}
\alias{skewProbit}

\title{Fitting Binary Regression with a Skew-Probit Link Function}

\description{
This function fits a binary regression with a skew-probit link function. 
Naive MLE, Jeffrey's prior and Cauchy prior type of penalization are implemented to find estimates.
}

\usage{
skewProbit(formula, data = list(), penalty = "Jeffrey", initial = NULL, 
	cvtCov = TRUE, delta0 = 3, level = 0.95)
}


\arguments{
  \item{formula}{an object of class "formula" as in \code{lm}, \code{glm}}
  \item{data}{an optional data frame, list or environment containing the variables in the model as in \code{lm}, \code{glm}}
  \item{penalty}{type of penalty function. Default option is "Jeffrey". 
	"Cauchy" will give estimates with Cauchy prior penaly function. "Naive" will give ML estimates.}
  \item{initial}{a logical value. If specified, it will be used for the initial value of numerical optimization.}
  \item{cvtCov}{a logical value. If it is true, then all numerical values will be standardized to have mean zero and unit standard deviation.}
  \item{delta0}{an initial guess of skewness parameter.}
  \item{level}{a confidence level. Default value is 0.95.}
}

\details{
This function uses \code{ucminf} package for optimization. Also package \code{sn} is necessary.
A detailed disscussion can be found in the reference below.
}
\value{
An object of class \code{skewProbit} is returned with 
\item{coefficients}{A named vector of coefficients}
\item{stderr}{Standard errors of coefficients}
\item{zscore}{Z-scores of coefficients}
\item{pval}{p-values of coefficients}
\item{lower}{Lower limits of confidence intervals}
\item{upper}{Upper limits of confidence intervals}
}
\references{
Identifiability and bias reduction in the skew-probit model for a binary response. 
\emph{To appear in Journal of Statistical Computation and Simulation}.
}
\author{
DongHyuk Lee, Samiran Sinha
}

\examples{

library(sn)
library(ucminf)

n <- 500
b0 <- 0.34
delta <- 4
b1 <- 1
b2 <- -0.7

set.seed(1234)
x1 <- runif(n, -2, 2)
x2 <- rnorm(n, sd = sqrt(4/3))
eta <- as.numeric(b0 + b1*x1 + b2*x2)
p <- psn(eta, alpha = delta)
y <- rbinom(n, 1, p)

\dontrun{
dat <- data.frame(y, x1 = x1, x2 = x2)
mod1 <- skewProbit(y ~ x1 + x2, data = dat, penalty = "Jeffrey", cvtCov = FALSE, level = 0.95)
mod2 <- skewProbit(y ~ x1 + x2, data = dat, penalty = "Naive", cvtCov = FALSE, level = 0.95)
mod3 <- skewProbit(y ~ x1 + x2, data = dat, penalty = "Cauchy", cvtCov = FALSE, level = 0.95)
summary(mod1)
summary(mod2)
summary(mod3)
}
}

\keyword{Skew-probit link }

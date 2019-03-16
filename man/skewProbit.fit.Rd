\name{skewProbit.fit}
\alias{skewProbit.fit}

\title{Fitting Binary Regression with a Skew-Probit Link Function}

\description{It is the default fitting method for \code{skewProbit}.}

\usage{
skewProbit.fit(y, x, penalty = "Jeffrey", initial = NULL, 
	cvtCov = TRUE, delta0 = 3, level = 0.95)
}

\arguments{
  \item{y}{a design matrix of dimension \code{n*p}}
  \item{x}{a vector of response of length \code{n}}
  \item{penalty}{type of penalty function. Default option is "Jeffrey". 
	"Cauchy" will give estimates with Cauchy prior penaly function. "Naive" will give ML estimates.}
  \item{initial}{a logical value. If specified, it will be used for the initial value of numerical optimization.}
  \item{cvtCov}{a logical value. If it is true, then all numerical values will be standardized to have mean zero and unit standard deviation.}
  \item{delta0}{an initial guess of skewness parameter.}
  \item{level}{a confidence level. Default value is 0.95.}
}

\value{
A list cotaining the following components:
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

x <- cbind(1, x1, x2)
\dontrun{
mod1 <- skewProbit.fit(y, x, penalty = "Jeffrey", cvtCov = FALSE)
mod2 <- skewProbit.fit(y, x, penalty = "Naive", cvtCov = FALSE)
mod3 <- skewProbit.fit(y, x, penalty = "Cauchy", cvtCov = FALSE)
mod1$coef
mod2$coef
mod3$coef

} 
}

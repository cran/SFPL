\name{sfpl}
\alias{sfpl}

\title{Sparse Fused Plackett-Luce}

\description{
Contains the main function of this package that is used to estimate the parameter of interest \eqn{\beta}. The inner workings of the function are described in Hermes et al., (2024).
}

\usage{
sfpl(x, y, ls_vec, lf_vec, epsilon, verbose)
}

\arguments{
  \item{x}{A \eqn{M \times p} matrix containing the values for the \eqn{p} objects variables across the \eqn{M} objects.}
  \item{y}{A list consisting of \eqn{K} matrices with each matrix containing (partial) rankings across \eqn{n} observations for group \eqn{k}.}
  \item{ls_vec}{Vector containing shrinkage parameters.}
  \item{lf_vec}{Vector containing fusion penalty parameters.}
  \item{epsilon}{Small positive value used to ensure that the penalty function is differentiable. Typically set at \eqn{10^{-5}}.}
  \item{verbose}{Boolean that returns the process of the parameter estimation.}
}

\value{
\item{beta_est}{A list of length ls_vec\eqn{\times}lf_vec that contains the parameter estimates \eqn{\hat{beta}} for each combination of ls_vec and lf_vec.}
}

\references{
1. Hermes, S., van Heerwaarden, J., and Behrouzi, P. (2024). Joint Learning from Heterogeneous Rank Data. arXiv preprint, arXiv:2407.10846 \cr
}

\author{
  Sjoerd Hermes\cr
  Maintainer: Sjoerd Hermes \email{sjoerd.hermes@wur.nl}
}

\examples{
\donttest{

# we first obtain the rankings and object variables
data(ghana)
y <- list(ghana[[1]], ghana[[2]])
x <- ghana[[3]]

# our next step consists of creating two vectors for the penalty parameters
ls_vec <- lf_vec <- c(0, 0.25)

# we choose epsilon to be small: 10^(-5), as we did in Hermes et al., (2024)
# now we can fit our model
epsilon <- 10^(-5)
verbose <- FALSE

result <- sfpl(x, y, ls_vec, lf_vec, epsilon, verbose)
}
}


\name{sfpl_select}
\alias{sfpl_select}

\title{Model selection for SFPL}

\description{
This function selects the "best" fitted SFPL model using either the AIC or the BIC, see Hermes et al., (2024).
}
\usage{
sfpl_select(beta_est, x, y, ls_vec, lf_vec)}
\arguments{
  \item{beta_est}{A list of length ls_vec\eqn{\times}lf_vec that contains the parameter estimates \eqn{\hat{beta}}, using either sfpl or sfpl_approx,
  for each combination of ls_vec and lf_vec.}
  \item{x}{A \eqn{M \times p} matrix containing the values for the \eqn{p} objects variables across the \eqn{M} objects.}
  \item{y}{A list consisting of \eqn{K} matrices with each matrix containing (partial) rankings across \eqn{n} observations for group \eqn{k}.}
  \item{ls_vec}{Vector containing shrinkage parameters.}
  \item{lf_vec}{Vector containing fusion penalty parameters.}
}

\value{
\item{model_aic}{A \eqn{p \times K} matrix containing the parameter estimates using the penalty parameters \eqn{\lambda_s, \lambda_f} as chosen by the AIC.}
\item{model_bic}{A \eqn{p \times K} matrix containing the parameter estimates using the penalty parameters \eqn{\lambda_s, \lambda_f} as chosen by the BIC.}
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

result <- sfpl_approx(x, y, ls_vec, lf_vec, epsilon, verbose)

# now we select the best models using our model selection function
sfpl_select(result, x, y, ls_vec, lf_vec)
}
}


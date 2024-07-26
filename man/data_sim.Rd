\name{data_sim}
\alias{data_sim}

\title{Rank data simulation}

\description{
Simulates (partial) rank data for multiple groups together with object variables.
}

\usage{
data_sim(m, M, n, p, K, delta, eta)
}

\arguments{
  \item{m}{Length of the partial ranking for each observation.}
  \item{M}{Total number of objects.}
  \item{n}{Number of observations (rankers) per group.}
  \item{p}{Number of object variables.}
  \item{K}{Number of groups.}
  \item{delta}{Approximate fraction of different coefficients across the \eqn{\beta^{(k)}}.}
  \item{eta}{Approximate fraction of sparse coefficients in \eqn{\beta^{(k)}} for all \eqn{k}.}
}

\value{
\item{y}{A list consisting of \eqn{K} matrices with each matrix containing (partial) rankings across \eqn{n} observations for group \eqn{k}.}
\item{x}{A \eqn{M \times p} matrix containing the values for the \eqn{p} objects variables across the \eqn{M} objects.}
\item{beta}{A \eqn{p \times K} matrix containing the true value of \eqn{\beta}, which was used to generate \eqn{y}.}
}

\references{
1. Hermes, S., van Heerwaarden, J., and Behrouzi, P. (2024). Joint Learning from Heterogeneous Rank Data. arXiv preprint, arXiv:2407.10846 \cr
}

\author{
  Sjoerd Hermes\cr
  Maintainer: Sjoerd Hermes \email{sjoerd.hermes@wur.nl}
}

\examples{
data_sim(3, 10, 50, 5, 2, 0.25, 0.25)
}

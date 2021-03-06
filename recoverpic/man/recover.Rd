% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{recover}
\alias{recover}
\title{Fill in the missing picture}
\usage{
recover(pic, lambda, ac = 1L)
}
\arguments{
\item{pic}{is the input matrix which needs to be full of}

\item{lambda}{is the penalty variable}

\item{ac}{is the bool variable, 1 means using accelerate method, 0 means no accelerate}
}
\value{
the matrix which has already full of
}
\description{
Fill in the missing picture
}
\examples{
require(recoverpic)
nimgac = recover(vimg,0.3)
}

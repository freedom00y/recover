\name{recoverpic}
\alias{recoverpic-package}
\alias{recoverpic}
\docType{package}
\title{
\packageTitle{recoverpic}
}
\description{
\packageDescription{recoverpic}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{recoverpic}
\packageIndices{recoverpic}
Using Accelerated Proximal Gradient method to fill in the missing value for picures

vimg = vague(img)

nimgac = recover(vimg,0.3)
}
\author{
\packageAuthor{recoverpic}

Maintainer: \packageMaintainer{recoverpic}
}
\references{
Beck, A., Teboulle, M. (2009) A Fast Iterative Shrinkage-Thresholding Algorithm for Linear Inverse Problems. SIAM J. IMAGING SCIENCES.
}

\examples{
vimg = vague(img)
nimgac = recover(vimg,0.3)
}

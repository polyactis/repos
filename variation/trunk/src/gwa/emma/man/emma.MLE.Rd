\name{emma.MLE}
\alias{emma.MLE}
\title{Maximum likelihood estimation using linear mixed model}
\description{
  Estimates maximum likelihood and its parameters using linear mixed
  model
}
\usage{
  emma.MLE (y, x, K, Z=NULL, ngrids=100, llim=-5, ulim=5,
  esp=1e-10, eig.L = NULL, eig.R = NULL)
}
\arguments{
  \item{y}{a size n vector containing response variables (or
    phenotypes), where n is the number of individuals}
  \item{x}{a size t vector containing indicator variables (or snps),
    where n is the number of individuals}
  \item{K}{A t by t matrix of kinship coefficients, representing the
    pairwise genetic relatedness between strains}
  \item{Z}{A n by t incidence matrix mapping each individual to a
    strain. If this is NULL, n and t should be equal and an identity
    matrix replace Z}
  \item{ngrids}{Number of grids to search optimal variance component}
  \item{llim}{Lower bound of log ratio of two variance components}
  \item{ulim}{Upper bound of log ratio of two variance components}
  \item{esp}{Tolerance of numerical precision error}
  \item{eig.L}{Eigenvector from K used in ML estimate. If
    specified, it may avoid redundant computation inside the function}  
  \item{eig.R}{Eigenvector from x and K used in ML estimate. If
    specified, it may avoid redundant computation inside the function}
}
\details{
  The following criteria must hold; otherwise an error occurs
  - [length of y] == [# rows in Z]
  - [length of x] == [# cols in Z] == [# rows in K] == [# cols in K]
  - rowSums(Z) should be a vector of ones
  - colSums(Z) should not contain zero elements
  - K must be a positive semidefinite matrix
}
\value{
  A list containing:
  \item{ML}{Maximum likleihood estimate of the data given the linear
    mixed model}
  \item{delta}{Ratio between genetic and random variance component}
  \item{vg}{The coefficients of genetic variance component}
  \item{ve}{The coefficients of random variance component}
}
\references{
  Kang HM, Zaitlen NA, Wade CM, Kirby A, Heckerman D, Daly MJ, and Eskin
  E, Efficient Control of Population Structure in Model Organism
  Association Mapping, Genetics 178:1709-1723, 2008
}
\author{Hyun Min Kang \email{h3kang@cs.ucsd.edu}}
\seealso{\code{\link{emma.ML.LRT}},\code{\link{emma.kinship}},\code{\link{emma.REMLE}}}
\examples{
  \dontrun{
    ## Load data
    data(emmadat)

    ## Run EMMA MLE
    rs <- emma.MLE(emmadat$ys[1,],emmadat$xs[1,],emmadat$K)

    ## return maximum likelihood
    rs$ML
  }
}
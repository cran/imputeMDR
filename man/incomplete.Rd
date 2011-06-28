\name{incomplete}
\alias{incomplete}
\docType{data}
\title{A simulated example data containing missing values}
\description{
    This data is an example of a simulated case-control study data with 10% missing genotype values.
    10 biallelic marker genotypes (ex : single nucleotide polymorphisms; SNPs) for 200 cases and 200 controls are included.
    Class variable (response; phenotype; disease status) is in the first column.
}
\usage{data(incomplete)}
\format{
  A data frame with 400 observations on the following 11 variables.
  \describe{
    \item{\code{class}}{a numeric vector}
    \item{\code{snp1}}{a numeric vector}
    \item{\code{snp2}}{a numeric vector}
    \item{\code{snp3}}{a numeric vector}
    \item{\code{snp4}}{a numeric vector}
    \item{\code{snp5}}{a numeric vector}
    \item{\code{snp6}}{a numeric vector}
    \item{\code{snp7}}{a numeric vector}
    \item{\code{snp8}}{a numeric vector}
    \item{\code{snp9}}{a numeric vector}
    \item{\code{snp10}}{a numeric vector}
  }
}

\references{
Namkung J, Elston RC, Yang JM, Park T. "Identification of gene-gene interactions in the presence of missing data using the multifactor dimensionality reduction method" Genet Epidemiol. 2009 Nov;33(7):646-56.
}
\examples{
data(incomplete)
}
\keyword{datasets}

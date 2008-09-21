\name{impute.mdr}
\alias{impute.mdr}
\alias{mdr}
\alias{MDR}
\alias{multifactor dimensionality reduction}
\title{The Multifactor Dimensionality Reduction (MDR) Analysis for Imcomplete Data}
\description{
This provides various approaches to handling missing values for the MDR analysis of incomplete data
}
\usage{
impute.mdr(dataset, colresp, cs, combi, cv.fold = 10, na.method = 0, randomize = FALSE)
}
\arguments{
  \item{dataset}{a matrix of SNP data with class variable. genotypes must be coded as allele counts (0,1,2).}
  \item{colresp}{ location of class variable in the dataset. no missing values are allowed for response variable}
  \item{cs}{ the value used to indicate "case" for class variable}
  \item{combi}{number of SNPs considered simultaneously as predictor variables. an order of interactions to analyze. }
  \item{cv.fold}{ number of fold in cross validation }
  \item{na.method}{ options for missing handling approaches. 
                       na.method = 0 for complete data, na.method = 1 for treating missing genotypes 
		       as another genotype category, na.method=2 for using available data for given 
		       number of SNPs under consideration as a model, na.method=3 for using method 
		       of imputing missing information by using EM (expectation-maximization) algorithm}
  \item{randomize}{logical. If 'TRUE' the cross validation sets are randomized }
}
\details{
}
\value{
  \item{min.comb }{combination with minimum error rate in each cross validation}
  \item{train.erate }{training error}
  \item{test.erate }{test error}
  \item{best.combi }{the best combination}
}
\references{}
\author{
Junghyun Namkung, Taeyoung Hwang, MinSeok Kwon, Sunggon Yi and Wonil Chung
Maintainer: Junghyun Namkung <jh.namkung@gmail.com>}
\note{ 
}
\examples{
## sample data with missing values
data(incomplete)
## analysis example of 2nd order gene-gene interaction test
impute.mdr(incomplete, colresp=1, cs=1, combi=2, cv.fold = 10,na.method=2)
}
\keyword{nonparametric}


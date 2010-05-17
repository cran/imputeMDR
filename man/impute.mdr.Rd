\name{impute.mdr}
\alias{impute.mdr}
\alias{mdr}
\alias{MDR}
\alias{multifactor dimensionality reduction}
\title{The Multifactor Dimensionality Reduction (MDR) Analysis for Imcomplete Data}
\description{
This function provides various approaches to handling missing values for the MDR analysis of incomplete data to identify gene-gene interactions using biallelic marker data in genetic association studies 
}
\usage{
impute.mdr(dataset, colresp, cs, combi, cv.fold = 10, na.method = 0, max_iter = 30, randomize = FALSE)
}
\arguments{
  \item{dataset}{A matrix of SNP data with class variable (response; phenotype; disease status). Genotypes must be coded as allele counts (0,1,2). Missing genotypes should be coded as 3}
  \item{colresp}{ Column number of class variables in the dataset. No missing value is allowed for the class variable}
  \item{cs}{ The value used to indicate "case (affected)" for class variable}
  \item{combi}{The number of SNPs considered simultaneously as predictor variables (An order of interactions to analyze) }
  \item{cv.fold}{ The number of folds k for k-fold cross-validation }
  \item{na.method}{ Options for missing handling approaches. 
                       na.method = 0 for complete data, na.method = 1 for treating missing genotypes 
		       as another genotype category, na.method=2 for using available data for given 
		       number of SNPs under consideration as a model, na.method=3 for using method 
		       of imputing missing information by using EM (expectation-maximization) algorithm}
  \item{max_iter}{The number of maximum iteration in EM impute approach (na.method=3). In order to apply one-step EM approach, set this argument as 1}
  \item{randomize}{Logical. If 'TRUE' the cross validation sets are randomized }
}
\details{
}
\value{
  \item{min.comb }{Marker combinations with the minimum error rate in each cross validation}
  \item{train.erate }{Training errors for selected marker combination}
  \item{test.erate }{Test error of the selected marker combination}
  \item{best.combi }{The best combination that was selected most frequently across k-fold cross validation}  
}
\references{Namkung J, Elston RC, Yang JM, Park T. "Identification of gene-gene interactions in the presence of missing data using the multifactor dimensionality reduction method" Genet Epidemiol. 2009 Nov;33(7):646-56.
 }
\author{
Junghyun Namkung, Taeyoung Hwang, MinSeok Kwon, Sunggon Yi and Wonil Chung

Maintainer: Junghyun Namkung <jh.namkung@gmail.com>
 }
\note{ 
}
\seealso{ }
\examples{
## sample data with missing values
data(incomplete)
## analysis example of 2nd order gene-gene interaction test
impute.mdr(incomplete, colresp=1, cs=1, combi=2, cv.fold = 10,na.method=2)
}
\keyword{nonparametric}


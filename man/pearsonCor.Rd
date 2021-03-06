% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bestKeeper.R
\name{pearsonCor}
\alias{pearsonCor}
\title{Analyzes pair-wise correlation}
\usage{
pearsonCor(expression, ctVal = TRUE)
}
\arguments{
\item{expression}{a matrix of expression levels. Each row corresponds to a sample and each column to a gene.}

\item{ctVal}{a logical value indicating data type. If ct-values are input, ctVal=TRUE, otherwise, ctVal=FALSE.}
}
\value{
A matrix of the Pearson correlation coefficient (r) and the probability p value.
}
\description{
This function performs numerous pair-wise correlation analyses of genes. 
Within each such correlation the Pearson correlation coefficient (r) and the probability p value are calculated.
}
\examples{
FIBct
pearsonCor(FIBct)
}
\references{
Pfaffl MW, Tichopad A, Prgomet C, Neuvians TP. Biotechnol Lett (2004) <doi: 10.1023/B:BILE.0000019559.84305.47>
}

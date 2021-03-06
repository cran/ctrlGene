% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bestKeeper.R
\name{cpSta}
\alias{cpSta}
\title{Calculates descriptive statistics}
\usage{
cpSta(expression, ctVal = TRUE)
}
\arguments{
\item{expression}{a matrix of expression levels. Each row corresponds to a sample and each column to a gene.}

\item{ctVal}{a logical value indicating data type. If ct-values are input, ctVal=TRUE, otherwise, ctVal=FALSE.}
}
\value{
A matrix of descriptive statistics:

N: number of samples;

GM[CP]: the geometric mean of CP;

AM[CP]: the arithmetic mean of CP;

Min[CP] and Max [CP]: the extreme values of CP;

SD[+/- CP]: the standard deviation of the CP;

CV[CP]: the coefficient of variance expressed as a percentage on the CP level;

Min[x-fold] and Max [x-fold]: the extreme values of expression levels expressed as an absolute x-fold over- or under-regulation coefficient;

SD[+/- x-fold]: standard deviation of the absolute regulation coefficients.
}
\description{
This function calculates descriptive statistics of genes.
}
\examples{
FIBct
cpSta(FIBct)
}
\references{
Pfaffl MW, Tichopad A, Prgomet C, Neuvians TP. Biotechnol Lett (2004) <doi: 10.1023/B:BILE.0000019559.84305.47>
}

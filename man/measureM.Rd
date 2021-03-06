% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geNorm.R
\name{measureM}
\alias{measureM}
\title{Calculates measure M}
\usage{
measureM(expression, ctVal = TRUE)
}
\arguments{
\item{expression}{a matrix of expression levels. Each row corresponds to a sample and each column to a gene.}

\item{ctVal}{a logical value indicating data type. If ct-values are input, ctVal=TRUE, otherwise, ctVal=FALSE.}
}
\value{
A sorted dataframe with two columns, 'Genes' and 'M' (the relative stability; lower means more stable).
}
\description{
This function calculates measure M according to algorithm of geNorm
}
\examples{
FIB
measureM(FIB,ctVal=FALSE)
FIBct
measureM(FIBct)
}
\references{
Vandesompele J, De Preter K, Pattyn F, Poppe B, Van Roy N, De Paepe A, Speleman F (2002) <doi: 10.1186/gb-2002-3-7-research0034>.
}

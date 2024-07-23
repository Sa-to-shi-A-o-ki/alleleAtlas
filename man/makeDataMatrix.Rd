\name{makeDataMatrix}
\alias{makeDataMatrix}
\title{makeDataMatrix}
\usage{
makeDataMatrix(data,alternative,alleleSequence,id.samples)
}
\description{
Make a list to be used for the model calculation.
}
\arguments{
\item{data}{A character vector loaded from a csv file.}
\item{alternative}{A numeric vector of size one, showing the number of allele types in the sample.}
\item{alleleSequence}{A character vector of allele sequences.}
\item{id.samples}{A numeric vector obtained by getIDsamples function.}
}
\value{
A list of numeric vectors.
}

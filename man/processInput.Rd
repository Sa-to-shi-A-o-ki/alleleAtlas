\name{processInput}
\alias{processInput}
\title{processInput}
\usage{
processInput(data)
}
\description{
Process input data.
}
\arguments{
\item{data}{A character vector read from a csv file. Each row in the csv file corresponds to each sampling location. The first and second columns are for the longitude and latitude, respectively. The third and more columns are for alleles.}
}
\value{
A list of the size four.
  The first element is a matrix of location data.
  The second element is a numeric vector showing the number of each allele type in the sample.
  The third element is a character vector of the allele sequences.
  The fourth element is a numeric vector of size 1. This shows the number of the allele types in the sample.
}

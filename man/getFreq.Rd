\name{getFreq}
\alias{getFreq}
\title{getFreq}
\usage{
getFreq(alternative,models,randomModel=NULL)
}
\description{
Calculate allele frequencies from the model objects.
}
\arguments{
\item{alternative}{A numeric vector of size one. The number of allele types.}
\item{models}{A list of calculated models.}
\item{randomModel}{A model used to calculate "models".}
}
\value{
A matrix of allele frequencies. The columns and rows correspond to the allele types and polygons, respectively.
}

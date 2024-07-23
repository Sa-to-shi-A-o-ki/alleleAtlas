\name{makeSampleFreq}
\alias{makeSampleFreq}
\title{makeSampleFreq}
\usage{
makeSampleFreq(dataMatrix,SPDF)
}
\description{
Calculate sample allele frequencies to be shown on a map.
}
\arguments{
\item{dataMatrix}{A list made by makeDataMatrix function.}
\item{SPDF}{A spatial polygons data frame object of meshes.}
}
\value{
A matrix of allele frequencies. The columns and rows correspond to the allele types and polygons, respectively.
}

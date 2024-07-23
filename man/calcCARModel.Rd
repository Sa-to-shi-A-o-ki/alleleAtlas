\name{calcCARModel}
\alias{calcCARModel}
\title{calcCARModel}
\usage{
calcCARModel(alleleSequence,randomModel,adjacency,dataMatrix,prior="jeffreystdf",verbose=FALSE)
}
\description{
Calculate conditional autoregression model to interpolate allele frequencies.
}
\arguments{
\item{alleleSequence}{A character vector of allele sequences.}
\item{randomModel}{Directly passed to "model" argument in f function.}
\item{adjacency}{A numeric matrix of adjacency matrix of polygons. Passed to "graph" argument in f function.}
\item{dataMatrix}{A list made by makeDataMatrix function.}
\item{prior}{Directly passed to prior in the "hyper" argument in f function.}
\item{verbose}{Directly passed to inla function. If TRUE, computation detail will be shown during calculation.}
}
\value{
A list of model calculation results. Each element of the list corresponds to the model of each allele type.
}

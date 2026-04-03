\name{calcCAR}
\alias{calcCAR}
\title{calcCAR}
\usage{
calcCAR(alternative,randomModel,adjacency,dataMatrix,method="group",prior="jeffreystdf",verbose=FALSE)
}
\description{
Calculate conditional autoregression model to interpolate allele frequencies.
}
\arguments{
\item{alternative}{An integer of the number of allele types.}
\item{randomModel}{Directly passed to "model" argument in f function.}
\item{adjacency}{A numeric matrix of adjacency matrix of polygons. Passed to "graph" argument in f function.}
\item{dataMatrix}{A list made by makeDataMatrix function.}
\item{method}{"group" or "replicate". Switches the method of estimation.}
\item{prior}{Directly passed to prior in the "hyper" argument in f function.}
\item{verbose}{Directly passed to inla function. If TRUE, computation detail will be shown during calculation.}
}
\value{
A list of model calculation results. Each element of the list corresponds to the model of each allele type.
}

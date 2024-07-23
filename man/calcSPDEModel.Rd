\name{calcSPDEModel}
\alias{calcSPDEModel}
\title{calcSPDEModel}
\usage{
calcSPDEModel(alleleSequence,alternative,dataMatrix,prior="flat",verbose=FALSE,mesh,data.loc)
}
\description{
Calculate SPDE model to interpolate allele frequencies.
}
\arguments{
\item{alleleSequence}{A character vector of allele sequences.}
\item{alternative}{An integer, indicating the number of allele types.}
\item{dataMatrix}{A list made by makeDataMatrix function.}
\item{prior}{Directly passed to prior in the "hyper" argument in f function.}
\item{verbose}{Directly passed to inla function. If TRUE, computation detail will be shown during calculation.}
\item{mesh}{A fm_mesh_2d object.}
\item{data.loc}{A matrix of the sampling locations, having two columns which corresponds to the longitude and latitude.}
}
\value{
A list of model calculation results. Each element of the list corresponds to the model of each allele type.
}

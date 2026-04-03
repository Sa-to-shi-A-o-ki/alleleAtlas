\name{calcSPDE}
\alias{calcSPDE}
\title{calcSPDE}
\usage{
calcSPDE(alternative,dataMatrix,mesh,data.loc,verbose=FALSE,method="group",spdeType="pc",sigma0=1,rangeDenom=5,range=c(NA,0.01),sigma=c(3,0.01))
}
\description{
Calculate SPDE model to interpolate allele frequencies.
}
\arguments{
\item{alternative}{An integer, indicating the number of allele types.}
\item{dataMatrix}{A list made by makeDataMatrix function.}
\item{mesh}{A fm_mesh_2d object.}
\item{data.loc}{A matrix of the estimating locations, having two columns which corresponds to the longitude and latitude.}
\item{verbose}{Directly passed to inla function. If TRUE, computation detail will be shown during calculation.}
\item{method}{"group" or "replicate". Switches the method of estimation.}
\item{spdeType}{"pc" or "matern". Switches the used model.}
\item{sigma0}{Used only in "matern" model. This is the marginal variance.}
\item{rangeDenom}{Used only in "matern" model. The range of estimating area / rangeDenom should be the distance at which the correlation of two points is ca. 0.13.}
\item{range}{Used only in "pc" model. Passed to prior.range in inla.spde2.pcmatern function.}
\item{sigma}{Used only in "pc" model. Passed to prior.sigma in inla.spde2.pcmatern function.}
}
\value{
A list of model calculation results. Each element of the list corresponds to the model of each allele type.
}

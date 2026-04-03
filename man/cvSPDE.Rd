\name{cvSPDE}
\alias{cvSPDE}
\title{cvSPDE}
\usage{
cvSPDE(data,SPDF,estimate.loc,mesh,allId.samples,uniqueAlleleIndex,sampleProb,method="group",spdeType="pc",sigma0=1,rangeDenom=5,range=c(NA,0.01),sigma=c(3,0.01))
}
\description{
Conduct leave one out cross validation for SPDE models.
}
\arguments{
\item{data}{A character vector loaded from a csv file.}
\item{SPDF}{A spatial polygons data frame object to calculate the frequencies.}
\item{estimate.loc}{A matrix of all the locations for estimation (not sampling location), having two columns which corresponds to the longitude and latitude.}
\item{mesh}{A fm_mesh_2d object.}
\item{allId.samples}{A numeric vector obtained by getIDsamples function for all locations.}
\item{uniqueAlleleIndex}{A list of location-specific alleles obtained by the fifth value of processInput function.}
\item{sampleProb}{A matrix of sample frequency. Can be obtained by makeSampleFreq function.}
\item{method}{Directly passed to method argument in calSPDE function.}
\item{spdeType}{Directly passed to spdeType argument in calcSPDE function.}
\item{sigma0}{Directly passed to sigma0 argument in calcSPDE function.}
\item{rangeDenom}{Directly passed to rangeDenom argument in calcSPDE function.}
\item{range}{Directly passed to range argument in calcSPDE function.}
\item{sigma}{Directly passed to sigma argument in calcSPDE function.}
}
\value{
A list of two numeric vectors. The first one is halved total absolute errors for each location. The second one is logloss for each location. Because logloss cannot be calculated on a location where location-specific allele is in it, logloss for such locations is always NA.
}

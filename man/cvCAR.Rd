\name{cvCAR}
\alias{cvCAR}
\title{cvCAR}
\usage{
cvCAR(data,randomModel,SPDF,data.loc,adjacency,allId.samples,uniqueAlleleIndex,sampleProb,method="group",prior="jeffreystdf")
}
\description{
Conduct leave one out cross validation for CAR models.
}
\arguments{
\item{data}{A character vector loaded from a csv file.}
\item{randomModel}{Directly passed to randomModel argument in calCARModel function.}
\item{SPDF}{A spatial polygons data frame object to calculate the frequencies.}
\item{data.loc}{A matrix of all the sampling locations, having two columns which corresponds to the longitude and latitude.}
\item{adjacency}{Directly passed to adjacency argument in calCARModel function.}
\item{allId.samples}{A numeric vector obtained by getIDsamples function for all locations.}
\item{uniqueAlleleIndex}{A list of location-specific alleles obtained by the fifth value of processInput function.}
\item{sampleProb}{A matrix of sample frequency. Can be obtained by makeSampleFreq function.}
\item{method}{Directly passed to method argument in calCARModel function.}
\item{prior}{Directly passed to prior argument in calcCARModel function.}
}
\value{
A list of two numeric vectors. The first one is halved total absolute errors for each location. The second one is logloss for each location. Because logloss cannot be calculated on a location where location-specific allele is in it, logloss for such locations is always NA.
}

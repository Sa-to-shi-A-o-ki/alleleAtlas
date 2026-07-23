\name{cvCAR}
\alias{cvCAR}
\title{cvCAR}
\usage{
cvCAR(data,randomModel,SPDF,adjacency,allId.samples,allDataMatrix,
allAlleleSequence,method="group",prior="jeffreystdf",infTreat=NA)
}
\description{
Conduct leave one out cross validation for CAR models.
}
\arguments{
\item{data}{A character vector loaded from a csv file.}
\item{randomModel}{Directly passed to randomModel argument in calCARModel function.}
\item{SPDF}{A spatial polygons data frame object to calculate the frequencies.}
\item{adjacency}{Directly passed to adjacency argument in calCARModel function.}
\item{allId.samples}{A numeric vector obtained by getIDsamples function for all locations.}
\item{allDataMatrix}{A list made by makeDataMatrix function for all locations.}
\item{allAlleleSequence}{A character vector of all allele sequences.}
\item{method}{Directly passed to method argument in calCARModel function.}
\item{prior}{Directly passed to prior argument in calcCARModel function.}
\item{infTreat}{How inf is treated in logloss and likelihood (mainly due to mesh-specific allele).}
}
\value{
A list of three numeric vectors.
The first one is total variation distance for each location.
The second one is logloss for each location. Because logloss cannot be calculated on a location where location-specific allele is in it, logloss for such locations is always NA.
The third one is negative log likelihood for each location. Just like logloss, it will be NA where location-specific allele is in the location.
}

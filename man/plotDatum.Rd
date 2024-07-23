\name{plotDatum}
\alias{plotDatum}
\title{plotDatum}
\usage{
plotDatum(polygon,data,title,markerLoc=NULL,tile="Esri.WorldGrayCanvas",col1="red",col2="yellow")
}
\description{
Plot an allele frequency or the other statistic as a heat map.
}
\arguments{
\item{polygon}{A spatial polygons data frame or fm_mesh_2d object used to calculate the frequencies.}
\item{data}{A numeric vector of data to be plotted.}
\item{title}{A character shown on the legend of the map.}
\item{markerLoc}{A numeric matrix of the sampling locations. It has two columns; the first is for the longitude, and the second for the latitude.}
\item{tile}{A character specifying the map. Please see provider argument in ?addProviderTiles.}
\item{col1}{A character specifying the heat map color for low values.}
\item{col2}{A character specifying the heat map color for high values.}
}
\value{
A leaflet map to be shown on a web browser.
}

\name{plotAllFreq}
\alias{plotAllFreq}
\title{plotAllFreq}
\usage{
plotAllFreq(polygon,data,groupList=NULL,dimx=10,dimy=10,type="pie",size=30,tile="Esri.WorldGrayCanvas",color=NULL)
}
\description{
Plot all frequencies on a map.
}
\arguments{
\item{polygon}{A spatial polygons data frame or fm_mesh_2d object used to calculate the frequencies.}
\item{data}{A matrix of allele frequencies.}
\item{groupList}{A list of integer vectors. Each element of the list constitutes an allele group which is shown in the same color on a map. The integers in the vector correspond to the index of allele types. So, for example, if there are 5 alleles and groupList = list(c(1),c(2,3),c(4,5)), 5 alleles are grouped into 3 allele groups. If one allele is included in multiple allele groups, the function returns an error. If one allele is not included in any allele groups, the function gives a warning but still works. When groupList=NULL (default option), no allele groups are made, and all allele types are shown on a map.}
\item{dimx}{An integer, specifying the number of columns in the plotting area. Only used when "polygon" is fm_mesh_2d object.}
\item{dimy}{An integer, specifying the number of rows in the plotting area. Only used when "polygon" is fm_mesh_2d object.}
\item{type}{A character, specifying the type of the charts. Please see ?addMinicharts.}
\item{size}{An integer, specifying the size of the charts.}
\item{tile}{A character specifying the map. Please see 'provider' argument in ?addProviderTiles.}
\item{color}{A vector of characters, specifying the color of the charts. When NULL (default), rainbow colors of 25\% transparency will be used.}
}
\value{
A leaflet map to be shown on a web browser.
}

\name{makeSpatialPolygons}
\alias{makeSpatialPolygons}
\title{makeSpatialPolygons}
\usage{
makeSpatialPolygons(autoWindow,lat,lon,minlon=-180,maxlon=180,minlat=-90,maxlat=90,data.loc=NULL,northPole=NULL,southPole=NULL,loop=NULL,forceMaxEast=NULL)
}
\description{
Make a spatial polygons data frame of rectangular meshes.
}
\arguments{
\item{autoWindow}{A logical vector of size one. Whether or not automatically define the area of meshes.}
\item{lat}{An interger. The number of vertical meshes.}
\item{lon}{An integer. The number of horizontal meshes.}
\item{minlon}{Necessary only when autoWindow=FALSE. An interger. The minimum longitude of area.}
\item{maxlon}{Necessary only when autoWindow=FALSE. An integer. The maximum longitude of area.}
\item{minlat}{Necessary only when autoWindow=FALSE. An interger. The minimum latitude of area.}
\item{maxlat}{Necessary only when autoWindow=FALSE. An integer. The maximum latitude of area.}
\item{data.loc}{Necessary only when autoWindow=TRUE. A matrix of coordinates used to define the area of meshes.}
\item{northPole}{A logical vector of size one. Whether or not all the meshes in the top row neighbor one another.}
\item{southPole}{A logical vector of size one. Whether or not all the meshes in the bottom row neighbor one another.}
\item{loop}{A boolean. Whether or not the leftmost and righrmost meshes neighbor one another.}
\item{forceMaxEast}{Necessary only when autoWindow=TRUE and the difference of the maximum and minimum longitudes is just 180 degree. A logical vector of size one. Use the eastmost (TRUE) or westmost  (FALSE) location as the rightmost window limit. In other words, when TRUE, 0 degree meridian is included in the area, and when FALSE, 180 degree meridian is included.}
}
\value{
A spatial polygons data frame.
}

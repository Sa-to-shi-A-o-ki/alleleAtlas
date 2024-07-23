\name{makeNeighbor}
\alias{makeNeighbor}
\title{makeNeighbor}
\usage{
makeNeighbor(path,lat,lon,northPole,southPole,loop)
}
\description{
Make a neighbor matrix of rectangular meshes and export it as a file.
}
\arguments{
\item{path}{A character. The path of exporting file.}
\item{lat}{An interger. The number of vertical meshes.}
\item{lon}{An integer. The number of horizontal meshes.}
\item{northPole}{A logical vector of size one. Whether or not all the meshes in the top row neighbor one another.}
\item{southPole}{A logical vector of size one. Whether or not all the meshes in the bottom row neighbor one another.}
\item{loop}{A logical vector of size one. Whether or not the leftmost and righrmost meshes neighbor one another.}
}
\value{
The function itself returns no value on R but exports a neighbor matrix as a file.
}

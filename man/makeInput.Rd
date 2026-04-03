\name{makeInput}
\alias{makeInput}
\title{makeInput}
\usage{
makeInput(lon,lat,fasta,corres)
}
\description{
Make a neighbor matrix of rectangular meshes and export it as a file.
}
\arguments{
\item{lon}{A vector of longitudes of sampling locations.}
\item{lat}{A vector of latitudes of sampling locations.}
\item{fasta}{A fasta file read by read.dna function in ape package.}
\item{corres}{A vector of characters, which can be used as an input of processInput function.}
}
\value{
The function itself returns no value on R but exports a neighbor matrix as a file.
}

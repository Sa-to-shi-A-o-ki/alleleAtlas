\name{getIDsamples}
\alias{getIDsamples}
\title{getIDsamples}
\usage{
getIDsamples(data.loc,SPDF,usest=FALSE)
}
\description{
Return a mesh id for each intersecting sampling location. Even when a sampling location intersects multiple meshes, this function returns only one id.
}
\arguments{
\item{data.loc}{A matrix of the sampling locations, having two columns which corresponds to the longitude and latitude.}
\item{SPDF}{A spatial polygons data frame object of meshes.}
\item{usest}{A logical vector of size one. If TRUE, use st_intersect function to obtain an intersection of the sampling location and meshes. If FALSE, use s2_intersection to obtain it. Try another if one does not work well. In some cases, both methods fails to generate proper IDs. In that case, please slightly change the meshes.}
}
\value{
A numeric vector indicating each sampling location intersects what mesh.
}

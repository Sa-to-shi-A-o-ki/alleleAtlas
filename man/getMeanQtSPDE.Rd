\name{getMeanQtSPDE}
\alias{getMeanQtSPDE}
\title{getMeanQtSPDE}
\usage{
getMeanQtSPDE(model,mesh,data.loc,n=1001,q=c(0.5,0.025,0.975))
}
\description{
Calculate mean and quantiles of allele frequencies from the SPDE model objects.
}
\arguments{
\item{model}{A list of calculated model.}
\item{mesh}{A mesh object used in the model.}
\item{data.loc}{A matrix of the estimating locations, having two columns which corresponds to the longitude and latitude.}
\item{n}{An integer of the resampling number. The larger number gives more precise result but takes longer time for calculation.}
\item{q}{A numeric vector to designate the quantiles. The default gives median and 95\% CI.}
}
\value{
A list of the allele frequency matrix of mean and quantiles. The columns and rows of the matrix correspond to the allele types and meshes, respectively.
}

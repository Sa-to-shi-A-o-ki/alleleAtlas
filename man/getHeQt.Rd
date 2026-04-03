\name{getHeQt}
\alias{getHeQt}
\title{getHeQt}
\usage{
getHeQt(model,n=1001,q=c(0.5,0.025,0.975))
}
\description{
Calculate quantiles of expected heterozygosity based on resamples from the posterior distribution.
}
\arguments{
\item{model}{A list of calculated models.}
\item{n}{An integer of the resampling number. The larger number gives more precise result but takes longer time for calculation.}
\item{q}{A numeric vector to designate the quantiles. The default gives median and 95\% CI.}
}
\value{
A matrix of excpected heterozygosity. The columns and rows of the matrix correspond to quantiles and meshes, respectively.
}

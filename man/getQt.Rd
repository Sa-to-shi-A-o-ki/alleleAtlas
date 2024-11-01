\name{getQt}
\alias{getQt}
\title{getQt}
\usage{
getQt(models,n=121,q=c(0.5,0.025,0.975))
}
\description{
Calculate quantiles of allele frequencies from the model objects.
}
\arguments{
\item{models}{A list of calculated models.}
\item{n}{An integer of the resampling number. The larger number gives more precise result but takes longer time for calculation.}
\item{q}{A numeric vector to designate the quantiles. The default gives median and 95\% CI.}
}
\value{
A list of the allele frequency matrix of each quantile. The columns and rows of the matrix correspond to the allele types and meshes, respectively.
}

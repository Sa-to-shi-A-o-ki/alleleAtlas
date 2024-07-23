\name{calcExpectedHeterozygosity}
\alias{calcExpectedHeterozygosity}
\title{calcExpectedHeterozygosity}
\usage{
calcExpectedHeterozygosity(freq)
}
\description{
Calculate expected heterozygosity from allele frequencies. The return value is not adjusted to the unbiased value by the sample size.
}
\arguments{
\item{freq}{A matrix of allele frequencies.}
}
\value{
A numeric vector of the expected heterozygosity for each location.
}

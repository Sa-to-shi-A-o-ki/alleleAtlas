\name{calcNucleotideDiversity}
\alias{calcNucleotideDiversity}
\title{calcNucleotideDiversity}
\usage{
calcNucleotideDiversity(freq,alleleSequence)
}
\description{
Calculate nucleotide diversity from allele frequencies. The return value is not adjusted to the unbiased value by the sample size.
}
\arguments{
\item{freq}{A matrix of allele frequencies.}
\item{alleleSequence}{A characteristic vector of nucleotide sequences of alleles. The nucleotide bases can be upper cases and lower cases.}
}
\value{
A numeric vector of the nucleotide diversity for each location.
}

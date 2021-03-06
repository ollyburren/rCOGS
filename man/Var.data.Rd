% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sCVPP.R
\name{Var.data}
\alias{Var.data}
\title{helper function to compute shrinkage factor for a quantitative trait}
\usage{
Var.data(f, N)
}
\arguments{
\item{f}{a vector of length $i$ of minor allele frequencies in controls matching the order of p}

\item{N}{a scalar or vector (of length $i$) of the total number of samples included in association test}
}
\value{
a vector of shrinkage factors
}
\description{
helper function to compute shrinkage factor for a quantitative trait
}

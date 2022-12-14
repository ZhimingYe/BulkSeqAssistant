% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Main.R
\name{MultiAOV}
\alias{MultiAOV}
\title{Multi group Parametric test}
\usage{
MultiAOV(Mat, FactorCol, Method = "LSD", OnlySig = F, pcutoff = 0.01)
}
\arguments{
\item{Mat}{Columns should be every variable to be tested, and another one column contains information about grouping, which is refer to the FactorCol}

\item{FactorCol}{Which column contains the grouping information. The number of it.}

\item{Method}{One of "LSD" or "SNK"}

\item{OnlySig}{Whether show only significant result.}

\item{pcutoff}{P value cut off.}
}
\description{
Multi group Parametric test
}
\author{
Zhiming Ye
}

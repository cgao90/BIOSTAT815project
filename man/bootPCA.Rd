% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{bootPCA}
\alias{bootPCA}
\title{An Rcpp function that carries out bootstrap PCA and calculats moment-based CIs and percentile CIs}
\usage{
bootPCA(Y, K, B, alpha, centerSamples = 0L)
}
\arguments{
\item{Y}{A (p x n) matrix of observed data points.}

\item{K}{num of leading PCs of interest}

\item{B}{number of bootstrap samples}

\item{alpha}{parameter for calculating percentile CIs}
}
\value{
A list of boostrap results
}
\description{
An Rcpp function that carries out bootstrap PCA and calculats moment-based CIs and percentile CIs
}

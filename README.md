# BIOSTAT815project
A reimplementation of bootSVD R package (Author: Dr.Aaron Fisher) in Rcpp

## Sample Code
### Install Packages 
```
install.packages("bootSVD")
install.packages("~/desktop/biostat815project_1.0.tar.gz", repos=NULL) 
```
### Load Packages
```
library(Rcpp)
library(RcppArmadillo)
library(bootSVD)
library(biostat815project)
```
### Generate simulated data
```
set.seed(0)
Y <- bootSVD::simEEG(n=100, centered=TRUE, wide=TRUE) 
Y_long = t(Y) # wide format to long format data
```

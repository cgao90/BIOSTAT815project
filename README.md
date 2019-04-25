# BIOSTAT815project
A reimplementation of bootSVD R package (Author: Dr.Aaron Fisher) in Rcpp
reference: (https://cran.r-project.org/web/packages/bootSVD/index.html)

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
### Implement bootstrap PCA
```
results <- biostat815project::bootPCA(Y_long,3,50,0.05)
```
### Plot results
```
k = 1
# bootstrapped PCs in low dimension indexed by K
AsByK <- results$full_LD_PC_dist
# bootstrapped PCs in high dimension indexed by K
VsByK = results$full_HD_PC_dist
# 95% CIs for PCs in low dimension
results$momentsAndMomentCI_LD$momentCI[[k]]
# 95% CIs for PCs in high dimension
results$momentsAndMomentCI_HD$momentCI[[k]]
# percentile CIs for PCs in low dimension
results$percentilesCI_LD[[k]]
# percentile CIs for PCs in high dimension
results$percentilesCI_HD[[k]]
```

# vistamilk

Predictions from an ensemble model: averages taken over predictions from:
* lasso from (glmnet library in R)
* random forest (ranger library in R) with and without regularization
* linear regression with 6 PCs with highest eigenvalues
* linear regression with 6 KPCs with highest eigenvalues using gaussian kernel
* random forest with 6 KPCs with highest eigenvalues using gaussian kernel
* random forest with a gaussian kernel
* support vector machine for regression (library e1071)
* pls with 3, 4, and 5 components
* linear regression with 15 waveles selected by a) clustering the wavelengths into 15 clusters and selecting the wavelength from each cluster that has the highest correlation with the response variable
* BART from bartMachine library

code at https://github.com/domijan/vistamilk


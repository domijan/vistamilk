# vistamilk competition

In this approach the training set was centered and scaled and the same transformation was applied to the test set. 
No outliers were removed and the full set of 1060 transmittance values was used for training and prediction. Casein micelle size was log transformed. Hierarchical clustering with $1 - $absolute value of correlation taken as distance matrix and Ward's clustering criterion was applied to the spectra. The resulting dendogram was visualised and the transmittance values were assigned to one of fifteen clusters.

The dataset was initially split into 50 random splits of training (50\%) and validation (50\%) sets to select the optimal algorithm for regression.

Following algorithms were used for regression:

* lasso from (glmnet library in R)
* random forest (ranger library in R) with and without regularization
* linear regression with 6 PCs with highest eigenvalues
* linear regression with 6 KPCs with highest eigenvalues using Gaussian kernel
* random forest with 6 KPCs with highest eigenvalues using Gaussian kernel
* random forest with a Gaussian kernel
* support vector machine for regression (library e1071)
* pls with 3, 4, and 5 components
* linear regression with 15 wavelengths selected by a) clustering the wavelengths into 15 clusters and selecting the wavelength from each cluster that has the highest correlation with the response variable
* BART from bartMachine library

The algorithms were tuned using further cross-validation of the training sets. 

Finally, an ensemble model was used: averages taken over predictions in test sets from the models described above.

The predictive ability of the algorithms was evaluated based on RMSE. The best scoring algorithms over the 50 random splits for all three traits was the ensamble model. 

The predictions of this model trained over the entire training set  were submitted to the competition. 



Files:

* analysis.R Code for the competition
* KDpresentation.Rmd: presentation for the workshop
* WorkshopPaperKD: participant section submission for the workshop paper
* WorkshopPaperKD.bib: bibliography

code at <https://github.com/domijan/vistamilk>


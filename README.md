# neighborhood-machine-learning

This repository contains code files to support the manuscript and supplementary files
"Variable selection in social-environmental data: Sparse regression and tree ensemble machine learning approaches"
Authors: Elizabeth Handorf, Yinuo Yin, Michael Slifker, and Shannon Lynch

All code is for binary outcomes, as code for linear models is very similar

Univariable regression, lasso: bfn_lso_bin_tt.R  
Sparse Group lasso: cluster_bin_tt.R  
Elastic Net: elnet_bin_tt.R  
Random Forests: RF_bin_tt.R   
Bagging: bagging_bin_tt.R  
Bayesian Adaptive Regression Trees: BART_bin_tt.R  

1.	Additional file 1
TIFF (.tiff) file
Title: Dendrogram for correlation between variables
Description: Dendrogram showing the relationships between the 1,000 elements of the covariate matrix X.  The horizontal red line represents a correlation of 0.8.
2.	Additional file 2
TIFF (.tiff) file
Title: Distribution of 10 variables associated with simulated outcomes
Description: Histograms showing the distributions of X1-X10, the elements of X used to simulate the outcomes Y.
3.	Additional file 3
Excel (.xlsx) file
Title: Detection rates for each variable
Description: Table showing the proportion of simulations where each varaible was selected (by each method)
4.	Additional file 4
TIFF (.tiff) file
Title: Correlation structure of 10 variables
Description: Correlations structure of X1-X10, the elements of X used to simulate the outcomes Y.  Blue represents a positive correlation and red a negative correlation, with darker colors indicating a stronger relationship.
5.	Additional file 5
Excel (.xlsx) file
Title: Full data results
Description: Table containing all finidings when the HCLST-CORR-SGL method was applied to the full PA PCa dataset (binary outcome: diagnosis with aggressive PCa)


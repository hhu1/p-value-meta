# p-value-meta
A script to combine p-value across multiple related traits. In a
typical genome-wide association study or a whole-exome sequencing
study, the investigator may want to test for associations between 
genes/SNVs and multiple related traits, e.g., smoking, blood pressure, 
BMI etc. This script allows you to combine the p-values across all tested 
traits under the null hypothesis that none of the trait is associated with 
the predictor (SNV or gene).

To use, load the function rvat_chi2 or brown into R and provide a p-value 
matrix as the argument. The rows of the matrix are the genes and the 
columns are the traits. The function will return a vector of p-values 
(one value for each gene).

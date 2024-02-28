# Tumor-Gene-Expression-Purity-Adjustment
Project: Purity adjusted subtyping scheme 

# Methods
1. Find expression data on Xena browser (Full expression data from CIT classifier + full expression matrix with probe intensities)
   
2. Pre-processing: Look at distribution, make PCA to test for batch effects, if any remove batch effect.
   
3. Subtype the breast cancer data based on already known data using classification models (kNN, distance to centroid, ssGSEA and potentially more)
   
4. ESTIMATE is used for predicting purity to determine if expression is driven by purity.
   
5. Adjust gene expression by purity by assuming 100% purity in all samples.
   
6. Run the CIT classifier (R package – CITBCMST)
   
7. Compare the classes from the articles training data and with purity adjusted data.
   
8. Discuss results
    

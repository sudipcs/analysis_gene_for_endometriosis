

https://doi.org/10.1016/B978-0-443-24028-7.00027-1


**1. endometriosisTrans.py**
Analyze endometriosis data to identify top marker genes using K-means, hierarchical clustering, and gene ranking techniques in Python

**2. LASSO_Coefficients.py** <br>

Feature Selection with LASSO: LASSO is applied with a dummy target (for demonstration) to select the most relevant genes. Here the genes with non-zero coefficients after LASSO are chosen. <br>
Clustering: The selected genes are used to cluster the samples. Then PCA is used for visualizing the clusters in 2D space. <br>
   
Adjusting the Regularization Parameter (Alpha): It may be necessary to adjust the regularization parameter to be less strict to ensure that LASSO retains more features. This can be achieved by either tuning the alpha value or manually selecting a lower value instead of relying on LassoCV's default cross-validation.

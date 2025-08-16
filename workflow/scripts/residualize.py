import numpy as np
from scipy import stats
from sklearn.linear_model import LinearRegression

def calculate_residual(phenotype_df, covariates_df, center=False):
    """
    Calculate normalized residual phenotypes using linear regression.
    
    This function replicates the TensorQTL Residualizer functionality:
    1. Fits linear regression: phenotype ~ covariates
    2. Returns residuals (phenotype - predicted_phenotype)
    3. Optionally centers and normalizes residuals
    
    Args:
        phenotype_df (pd.DataFrame): Expression data (genes x samples)
        covariates_df (pd.DataFrame): Covariates data (samples x covariates)
        center (bool): Whether to center and normalize residuals
    
    Returns:
        np.ndarray: Residualized phenotypes (genes x samples)
    """
    # Convert to numpy arrays
    phenotype_values = phenotype_df.values  # genes x samples
    covariates_values = covariates_df.values  # samples x covariates
    
    # Ensure covariates has intercept (add column of 1s if not present)
    if not np.allclose(covariates_values[:, 0], 1.0):
        covariates_with_intercept = np.column_stack([np.ones(covariates_values.shape[0]), covariates_values])
    else:
        covariates_with_intercept = covariates_values
    
    # Initialize output array
    residuals = np.zeros_like(phenotype_values)
    
    # For each gene, fit linear regression and get residuals
    for i in range(phenotype_values.shape[0]):
        # Fit linear regression: phenotype[i] ~ covariates
        lr = LinearRegression(fit_intercept=False)  # intercept already in covariates
        lr.fit(covariates_with_intercept, phenotype_values[i, :])
        
        # Get residuals: actual - predicted
        predicted = lr.predict(covariates_with_intercept)
        residuals[i, :] = phenotype_values[i, :] - predicted
    
    # Center and normalize if requested
    if center:
        # Center: subtract mean across samples for each gene
        residuals = residuals - np.mean(residuals, axis=1, keepdims=True)
        
        # Normalize: divide by standard deviation across samples for each gene
        std_dev = np.std(residuals, axis=1, keepdims=True)
        # Avoid division by zero
        std_dev[std_dev == 0] = 1.0
        residuals = residuals / std_dev
    
    return residuals
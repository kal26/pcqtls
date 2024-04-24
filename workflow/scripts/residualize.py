import tensorqtl
import torch
from tensorqtl.core import Residualizer

# adapted from tensorqtl
def calculate_residual(phenotype_df, covariates_df, center=False):
    """Calculate normalized residual phenotypes"""
    # set up residualize device
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    residualizer = Residualizer(torch.tensor(covariates_df.values, dtype=torch.float32).to(device))
    phenotype_t = torch.tensor(phenotype_df.values, dtype=torch.float).to(device)

    # residualize
    phenotype_res_t = residualizer.transform(phenotype_t)  # phenotypes x samples

    # center and normalize
    if center:
        phenotype_res_t = tensorqtl.core.center_normalize(phenotype_res_t, dim=1)

    return phenotype_res_t.numpy()
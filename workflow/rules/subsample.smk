"""
Subsampling Rules

This module contains rules for creating subsamples of expression and covariate data
for testing and validation purposes.
"""

rule get_subsample:
    """
    Create subsamples of expression and covariate data for testing.
    
    This rule generates smaller datasets by selecting a subset of samples from the
    full expression and covariate data, either using random sampling or taking
    the first N samples.
    """
    input: 
        expression = f"{config['expression_dir']}{{TISSUE}}.v8.normalized_expression.bed",
        covariates = f"{config['covariates_dir']}{{TISSUE}}.v8.covariates.txt"
    
    output:
        subsample_expression = f"output/subsample/normalized_expression/{{TISSUE}}.v8.normalized_expression.bed",
        subsample_covariates = f"output/subsample/covariates/{{TISSUE}}.v8.covariates.txt"
    
    params:
        use_scramble_order = False,
        num_samples = None,
        code_dir = config['code_dir']
    
    resources:
        mem = "8G",
        time = "0:30:00"
    
    threads: 1
    
    shell:
        """
        python {params.code_dir}/make_subsamples.py \
            --expression {input.expression} \
            --covariates {input.covariates} \
            --output-expression {output.subsample_expression} \
            --output-covariates {output.subsample_covariates} \
            --use-scramble
        """
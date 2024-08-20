rule get_subsample:
    input: 
        normalized_expression = expression_dir + '{TISSUE}.v8.normalized_expression.bed',
        full_covariates = covariates_dir + '{TISSUE}.v8.covariates.txt'
    params:
        use_scramble_order = False
    conda:
        'tensorqtl_r'
    output:
        subsample_expression = subsample_dir + 'normalized_expression/{TISSUE}.v8.normalized_expression.bed',
        subsample_covarience = subsample_dir + 'covariates/{TISSUE}.v8.covariates.txt
    script:
        '../scripts/make_subsamples.py'
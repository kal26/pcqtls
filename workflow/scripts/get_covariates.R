#------------------------------------------------------------------#
# Generate a file with conditions based on user-defined:           #
#     covars_yn - 'yes_covars' OR 'no_covars' OR 'both_covars'     #
#     pc_or_peer - 'pcs_only' OR 'peers_only'                      #
#     n_factors_max - [integer]                                    #
#     factors_break - [integer]                                    #
#     output_file - [character]                                    #
#------------------------------------------------------------------#

library(dplyr)
library(stringr)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
  
} else if (length(args) > 1) {
  
  covars_yn = as.character(args[1])
  pc_or_peer = as.character(args[2])
  n_factors_max = as.numeric(args[3])
  factors_break = as.numeric(args[4])
  output_file_path = as.character(args[5])
  
}


#------------#
# FUNCTIONS  #
#------------#
make_conditions_df = function(covars_yn, pc_or_peer, n_factors_max, factors_break){
  
  if (covars_yn %in% c('yes_covars', 'no_covars')){
    n_covar_conditions = 1 
  } else if (covars_yn == 'both_covars'){
    n_covar_conditions = 2
  }
    
  # IF only PCs
  if (pc_or_peer == 'pcs_only'){
      
    pcs = seq(0, n_factors_max, by=factors_break)
      
    conditions_df = data.frame(matrix(nrow=n_covar_conditions*length(pcs), ncol=3))
    names(conditions_df) = c('Covariates', 'N_PCs', 'N_PEERs')
      
    conditions_df$N_PCs = pcs
    conditions_df$N_PEERs = 0
    
    if (covars_yn %in% c('yes_covars', 'no_covars')){
      conditions_df$Covariates = covars_yn
    } else if (covars_yn == 'both_covars'){
      conditions_df$Covariates = c(rep('yes_covars', nrow(conditions_df)/2),
                                   rep('no_covars', nrow(conditions_df)/2))
    }
      
  # IF only PEER factors
  } else if (pc_or_peer == 'peers_only'){
      
    peers = seq(0, n_factors_max, by=factors_break)
      
    conditions_df = data.frame(matrix(nrow=n_covar_conditions*length(peers), ncol=3))
    names(conditions_df) = c('Covariates', 'N_PCs', 'N_PEERs')
      
    conditions_df$Covariates = 'yes'
    conditions_df$N_PCs = 0
    conditions_df$N_PEERs = peers
    
    if (covars_yn %in% c('yes_covars', 'no_covars')){
      conditions_df$Covariates = covars_yn
    } else if (covars_yn == 'both_covars'){
      conditions_df$Covariates = c(rep('yes_covars', nrow(conditions_df)/2),
                                   rep('no_covars', nrow(conditions_df)/2))
    }
      
  }
    
  return(conditions_df)
  
}

#------------------------------#
# MAKE CONDITIONS & SAVE FILE  #
#------------------------------#
conditions_df = make_conditions_df(covars_yn, pc_or_peer, n_factors_max, factors_break)

write.csv(conditions_df, output_file_path,
          row.names=FALSE)

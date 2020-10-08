# calculate variance of dirichlet multinomial for one sample for one cluster
var_dirichlet_multinomial <- function(alphas,ind,size){
  sum_alphas <- sum(alphas)
  prop_alphas <- alphas[ind]/sum_alphas
  return(size*prop_alphas*(1-prop_alphas)*((size+sum_alphas)/(1+sum_alphas)))
}
# calculate covariance of DM for one sample between two clusters
covar_dirichlet_multinomial <- function(alphas,ind1,ind2,size){
  sum_alphas <- sum(alphas)
  prop_alphas_1 <- alphas[ind1]/sum_alphas
  prop_alphas_2 <- alphas[ind2]/sum_alphas
  return(size*prop_alphas_1*prop_alphas_2*((size+sum_alphas)/(1+sum_alphas)))
}

# covariance matrix for one sample
cov_matrix_dirichlet_multinomial <- function(alphas,size){
  cov_mat <- matrix(NA,nrow = length(alphas),ncol=length(alphas),dimnames = list(paste0("a_",round(alphas,2)),paste0("a_",round(alphas,2))))
  bool_ut <- upper.tri(cov_mat)
  for( j in seq_len(length(alphas)-1)){
    cov_mat[j,bool_ut[j,]] <- sapply((j+1):length(alphas),function(i) covar_dirichlet_multinomial(alphas,i,j,size))
  }
  cov_mat[lower.tri(cov_mat)] <- t(cov_mat)[lower.tri(cov_mat)]
  for( j in seq_along(alphas)){
    cov_mat[j,j] <- var_dirichlet_multinomial(alphas,j,size)
  }
  cov_mat
}

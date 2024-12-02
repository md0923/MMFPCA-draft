###########################################################################################
## Description: Functions for generating multilevel multivariate functional data sets with 
##              shared or correlated eigenscores described in 'Joint Modeling of Evoked and 
##              Induced Event-Related Spectral Perturbations.'
###########################################################################################
## Functions included:
## Main function:
##    1. shared_gen: Function that generates a multilevel multivariate functional data set
##                   with shared eigenscores across different variates
##    2. corr_gen: Function that generates a multilevel multivariate functional data set
##                 with eigenscores correlated between different variates
###########################################################################################


shared_gen = function(
  eigen_lvl1_var1,                  # subject-level eigenfunctions for variate 1 (matrix of dimension M*c1)
  eigen_lvl1_var2,                  # subject-level eigenfunctions for variate 2 (matrix of dimension M*c1)
  eigen_lvl2_var1,                  # trial-level eigenfunctions for variate 1 (matrix of dimension M*c2)
  eigen_lvl2_var2,                  # trial-level eigenfunctions for variate 2 (matrix of dimension M*c2)
  z_mu_var1,                        # overall mean function for variate 1 (array of length M)
  z_mu_var2,                        # overall mean function for variate 2 (array of length M)
  U,                                # subject-level eigenvalues (array of length c1)
  V,                                # trial-level eigenvalues (array of length c2)
  sigma2_e,                         # variance of measurement error (scalar)
  N,                                # number of subjects (scalar)
  R){                               # number of trials per subject (scalar)
  
  #########################################################################################
  ## Description: Function constructing the multilevel multi- and uni-variate two-dimensional 
  ##              eigenfunctions described in the simulation section of the paper.
  ## Definition   M: number of sampling points (same for all variates)
  ##              c1: number of subject-level multivariate eigenfunctions 
  ##              c2: number of trial-level multivariate eigenfunctions
  ##              N: number of subjects
  ##              R: number of trials per subject
  ## Args:        see above
  ## Returns:     a list of generated multilevel multivariate functional data
  ##              z_var1: multilevel functional data for variate 1 (matrix of dimension NR*M)
  ##              z_var2: multilevel functional data for variate 2 (matrix of dimension NR*M)
  ##              id_array: id information for each trial to identify the subject it belong to (array of length NR)
  #########################################################################################
  
  M = nrow(eigen_lvl1_var1)         # number of sampling grids
  c1 = length(U)                    # number of subject-level eigencomponents
  c2 = length(V)                    # number of trial-level eigencomponents
  
  # independently sample shared subject- and trial-level eigenscores from Gaussian distributions
  subject_score = rmvnorm(N, sigma = diag(U))
  trial_score = rmvnorm(N*R, sigma = diag(V))
  
  # construct multilevel multivariate functional data using K-L decomposition
  z_var1 = z_var2 = matrix(nrow = N*R, ncol = M)
  id_array = array(dim = N*R)
  for (i in 1:N){
    for (r in 1:R){
      index = R*(i-1) + r
      # variate 1
      z_var1[index,] = z_mu_var1 +                            # overall mean function 
        subject_score[i,] %*% t(eigen_lvl1_var1) +            # subject-specific deviation from mean
        trial_score[index,] %*% t(eigen_lvl2_var1) +          # trial-specific deviation from subject-mean
        rnorm(M, sd = sqrt(sigma2_e))                         # independent measurement error
      # variate 2
      z_var2[index,] = z_mu_var2 +                            # overall mean function  
        subject_score[i,] %*% t(eigen_lvl1_var2) +            # subject-specific deviation from mean
        trial_score[index,] %*% t(eigen_lvl2_var2) +          # trial-specific deviation from subject-mean
        rnorm(M, sd = sqrt(sigma2_e))                         # independent measurement error
      # save the subject information for the generated trial
      id_array[index] = i
    }
  }
  
  # assemble generated multilevel multivariate functional data and id array into one list
  z_data = list(z_var1 = z_var1,
                z_var2 = z_var2,
                id_array = id_array)
  return(z_data)
}




corr_gen = function(
  eigen_lvl1_var1,                  # subject-level eigenfunctions for variate 1 (matrix of dimension M*c1)
  eigen_lvl1_var2,                  # subject-level eigenfunctions for variate 2 (matrix of dimension M*c1)
  eigen_lvl2_var1,                  # trial-level eigenfunctions for variate 1 (matrix of dimension M*c2)
  eigen_lvl2_var2,                  # trial-level eigenfunctions for variate 2 (matrix of dimension M*c2)
  z_mu_var1,                        # overall mean function for variate 1 (array of length M)
  z_mu_var2,                        # overall mean function for variate 2 (array of length M)
  U,                                # subject-level eigenvalues (array of length c1)
  V,                                # trial-level eigenvalues (array of length c2)
  rho_lvl1,                         # correlation between subject-level scores (scalar < 1)
  rho_lvl2,                         # correlation between trial-level scores (scalar < 1)
  sigma2_e,                         # variance of measurement error (scalar)
  N,                                # number of subjects (scalar)
  R){                               # number of trials per subject (scalar)
  
  #########################################################################################
  ## Description: Function constructing the multilevel multi- and uni-variate two-dimensional 
  ##              eigenfunctions described in the simulation section of the paper.
  ## Definition   M: number of sampling points (same for all variates)
  ##              c1: number of subject-level multivariate eigenfunctions 
  ##              c2: number of trial-level multivariate eigenfunctions
  ##              N: number of subjects
  ##              R: number of trials per subject
  ## Args:        see above
  ## Returns:     a list of generated multilevel multivariate functional data
  ##              z_var1: multilevel functional data for variate 1 (matrix of dimension NR*M)
  ##              z_var2: multilevel functional data for variate 2 (matrix of dimension NR*M)
  ##              id_array: id information for each trial to identify the subject it belong to (array of length NR)
  #########################################################################################
  
  M = nrow(eigen_lvl1_var1)         # number of sampling grids
  c1 = length(U)                    # number of subject-level eigencomponents
  c2 = length(V)                    # number of trial-level eigencomponents
  
  # specify subject and trial-level eigenscore covariance matrices (block matrices)
  joint_cov_subject = diag(rep(U,2))
  joint_cov_subject[(c1+1):(2*c1), 1:c1] = joint_cov_subject[1:c1, (c1+1):(2*c1)] = rho_lvl1 * diag(U)
  joint_cov_trial = diag(rep(V,2))
  joint_cov_trial[(c2+1):(2*c2), 1:c2] = joint_cov_trial[1:c2, (c2+1):(2*c2)] = rho_lvl2 * diag(V)
  
  # generate subject- and trial-level eigenscores that are correlated between two variates
  subject_score = rmvnorm(N, sigma = joint_cov_subject)
  trial_score = rmvnorm(N*R, sigma = joint_cov_trial)
  subject_score_var1 = subject_score[,1:c1]
  subject_score_var2 = subject_score[,(c1+1):(2*c1)]
  trial_score_var1 = trial_score[,1:c2]
  trial_score_var2 = trial_score[,(c2+1):(2*c2)]
  
  # construct multilevel multivariate functional data using K-L decomposition
  z_var1 = z_var2 = matrix(nrow = N*R, ncol = M)
  id_array = array(dim = N*R)
  for (i in 1:N){
    for (r in 1:R){
      index = R*(i-1) + r
      # variate 1
      z_var1[index,] = z_mu_var1 +                            # overall mean function 
        subject_score[i,] %*% t(eigen_lvl1_var1) +            # subject-specific deviation from mean
        trial_score[index,] %*% t(eigen_lvl2_var1) +          # trial-specific deviation from subject-mean
        rnorm(M, sd = sqrt(sigma2_e))                         # independent measurement error
      # variate 2
      z_var2[index,] = z_mu_var2 +                            # overall mean function  
        subject_score[i,] %*% t(eigen_lvl1_var2) +            # subject-specific deviation from mean
        trial_score[index,] %*% t(eigen_lvl2_var2) +          # trial-specific deviation from subject-mean
        rnorm(M, sd = sqrt(sigma2_e))                         # independent measurement error
      # save the subject information for the generated trial
      id_array[index] = i
    }
  }
  
  # assemble generated multilevel multivariate functional data and id array into one list
  z_data = list(z_var1 = z_var1,
                z_var2 = z_var2,
                id_array = id_array)
  return(z_data)
}

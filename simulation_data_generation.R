shared_gen = function(
  eigen_lvl1_var1, # Txc1: each column is a subject-level eigenfunction for variate 1
  eigen_lvl1_var2, # Txc1: each column is a subject-level eigenfunction for variate 2 
  eigen_lvl2_var1, # Txc2: each column is a trial-level eigenfunction for variate 1 
  eigen_lvl2_var2, # Txc2: each column is a trial-level eigenfunction for variate 2
  z_mu_var1, # Tx1: overall mean for variate 1
  z_mu_var2, # Tx1: overall mean for variate 2
  U, # c1x1: subject-level eigenvalues
  V, # c2x1: trial-level eigenvalues
  sigma2_e, # scalar: variance of measurement error
  N, # scalar: number of subjects
  R){# scalar: number of trials per subject
  n_obs_1 = nrow(eigen_lvl1_var1)
  n_obs_2 = nrow(eigen_lvl1_var2)
  c1 = length(U)
  c2 = length(V)
  # generate shared eigenscores
  subject_score = rmvnorm(N, sigma = diag(U))
  trial_score = rmvnorm(N*R, sigma = diag(V))
  # generate simulated data with measurement errors
  z_var1 = matrix(nrow = N*R, ncol = n_obs_1)
  z_var2 = matrix(nrow = N*R, ncol = n_obs_2)
  for (i in 1:N){
    for (r in 1:R){
      index = R*(i-1) + r
      z_var1[index,] = z_mu_var1 + subject_score[i,] %*% t(eigen_lvl1_var1) +
        trial_score[index,] %*% t(eigen_lvl2_var1) +
        rnorm(n_obs_1, sd = sqrt(sigma2_e))
      z_var2[index,] = z_mu_var2 + subject_score[i,] %*% t(eigen_lvl1_var2) +
        trial_score[index,] %*% t(eigen_lvl2_var2) +
        rnorm(n_obs_2, sd = sqrt(sigma2_e))
    }
  }
  z_data = list(z_var1 = z_var1,
                z_var2 = z_var2)
  return(z_data)
}



corr_gen = function(
  eigen_lvl1_var1, # Txc1: each column is a subject-level eigenfunction for variate 1
  eigen_lvl1_var2, # Txc1: each column is a subject-level eigenfunction for variate 2 
  eigen_lvl2_var1, # Txc2: each column is a trial-level eigenfunction for variate 1 
  eigen_lvl2_var2, # Txc2: each column is a trial-level eigenfunction for variate 2
  z_mu_var1, # Tx1: overall mean for varaite 1
  z_mu_var2, # Tx1: overall mean for variate 2
  U, # c1x1: subject-level eigenvalues
  V, # c2x1: trial-level eigenvalues
  rho_lvl1, # scalar: the correlation between subject-level scores of the two variates
  rho_lvl2, # scalar: the correlation between trial-level scores of the two variates
  sigma2_e, # scalar: variance of measurement error
  N, # scalar: number of subjects
  R){ # percentage of variance explained threshold for prediction
  n_obs_1 = nrow(eigen_lvl1_var1)
  n_obs_2 = nrow(eigen_lvl1_var2)
  c1 = length(U)
  c2 = length(V)
  # generate subject and trial-level score covariance matrices that make var1 and var2 scores correlated
  require(mvtnorm)
  joint_cov_subject = diag(rep(U,2))
  joint_cov_subject[(c1+1):(2*c1), 1:c1] = joint_cov_subject[1:c1, (c1+1):(2*c1)] = rho_lvl1 * diag(U)
  joint_cov_trial = diag(rep(V,2))
  joint_cov_trial[(c2+1):(2*c2), 1:c2] = joint_cov_trial[1:c2, (c2+1):(2*c2)] = rho_lvl2 * diag(V)
  # generate individual eigenscores that are correlated between two variates
  subject_score = rmvnorm(N, sigma = joint_cov_subject)
  trial_score = rmvnorm(N*R, sigma = joint_cov_trial)
  subject_score_var1 = subject_score[,1:c1]
  subject_score_var2 = subject_score[,(c1+1):(2*c1)]
  trial_score_var1 = trial_score[,1:c2]
  trial_score_var2 = trial_score[,(c2+1):(2*c2)]
  # generate simulated data with measurement errors
  z_var1 = matrix(nrow = N*R, ncol = n_obs_1)
  z_var2 = matrix(nrow = N*R, ncol = n_obs_2)
  for (i in 1:N){
    for (r in 1:R){
      index = R*(i-1) + r
      z_var1[index,] = z_mu_var1 + subject_score_var1[i,] %*% t(eigen_lvl1_var1) +
        trial_score_var1[index,] %*% t(eigen_lvl2_var1) +
        rnorm(n_obs_1, sd = sqrt(sigma2_e))
      z_var2[index,] = z_mu_var2 + subject_score_var2[i,] %*% t(eigen_lvl1_var2) +
        trial_score_var2[index,] %*% t(eigen_lvl2_var2) +
        rnorm(n_obs_2, sd = sqrt(sigma2_e))
    }
  }
  z_data = list(z_var1 = z_var1,
                z_var2 = z_var2)
  return(z_data)
}

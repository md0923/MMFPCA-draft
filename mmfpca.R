mmfpca = function(
    z_var1,
    z_var2,
    id_var1,
    id_var2,
    x1, y1, x2, y2, # the axises of two dimensionals within each variate
    mufpca_pve = 0.99, # percentage of variance explained threshold for fitting 
    pred_pve){ # percentage of variance explained threshold for prediction
  # estimate the overall mean
  z_var1_mean = apply(z_var1, 2, mean)
  var1_mean_dt = data.frame(x1 = rep(x1, each = length(y1)),
                            y1 = rep(y1, length(x1)),
                            z_mean = c(z_var1_mean))
  var1_mu_fit = as.vector(predict(gam(z_mean ~ s(x1, y1), data = var1_mean_dt), 
                                  newdata = data.frame(x1 = rep(x1, each = length(y1)),
                                                       y1 = rep(y1, length(x1)))))
  z_var2_mean = apply(z_var2, 2, mean)
  var2_mean_dt = data.frame(x2 = rep(x2, each = length(y2)),
                            y2 = rep(y2, length(x2)),
                            z_mean = c(z_var2_mean))
  var2_mu_fit = as.vector(predict(gam(z_mean ~ s(x2, y2), data = var2_mean_dt), 
                                  newdata = data.frame(x2 = rep(x2, each = length(y2)),
                                                       y2 = rep(y2, length(x2)))))
  z_uni_pred_var1 = z_multi_pred_var1 = z_var1
  z_uni_pred_var2 = z_multi_pred_var2 = z_var2
  # fit multilevel univariate FPCA for each single variates
  mod_var1 = mfpca.face_center(Y = z_var1,
                               id = rep(1:N, each = R),
                               argvals = c(1:n_obs_1),
                               twoway = FALSE,
                               p = 4,
                               pve = mufpca_pve,
                               knots = floor(n_obs_1/5),
                               mu = var1_mu_fit)
  mod_var2 = mfpca.face_center(Y = z_var2,
                               id = rep(1:N, each = R),
                               argvals = c(1:n_obs_2),
                               twoway = FALSE,
                               p = 4,
                               pve = mufpca_pve,
                               knots = floor(n_obs_2/5),
                               mu = var2_mu_fit)
  # calculate the multivariate eigen components using univariate eigenscores
  m1_lvl1 = mod_var1$npc$level1
  m1_lvl2 = mod_var1$npc$level2
  m2_lvl1 = mod_var2$npc$level1
  m2_lvl2 = mod_var2$npc$level2
  m_sum_lvl1 = m1_lvl1 + m2_lvl1
  m_sum_lvl2 = m1_lvl2 + m2_lvl2
  # subject-level multivariate eigen components
  joint_score_lvl1 = cbind(mod_var1$scores$level1[,1:m1_lvl1],
                           mod_var2$scores$level1[,1:m2_lvl1])
  Z_joint_lvl1 = t(joint_score_lvl1) %*% joint_score_lvl1 / (N-1)
  Z_joint_lvl1.eigen = eigen(Z_joint_lvl1)
  evalue_lvl1_multi_est = Z_joint_lvl1.eigen$values
  eigen_lvl1_var1_multi_est = mod_var1$efunctions$level1[,1:m1_lvl1] %*%
    Z_joint_lvl1.eigen$vectors[1:m1_lvl1, ]
  eigen_lvl1_var2_multi_est = mod_var2$efunctions$level1[,1:m2_lvl1] %*%
    Z_joint_lvl1.eigen$vectors[c(m1_lvl1 + 1:m2_lvl1), ]
  score_lvl1_multi_est = mod_var1$scores$level1[,1:m1_lvl1] %*%
    Z_joint_lvl1.eigen$vectors[1:m1_lvl1, ] +
    mod_var2$scores$level1[,1:m2_lvl1] %*%
    Z_joint_lvl1.eigen$vectors[c(m1_lvl1 + 1:m2_lvl1), ]
  # trial-level multivariate eigen components
  joint_score_lvl2 = cbind(mod_var1$scores$level2[,1:m1_lvl2],
                           mod_var2$scores$level2[,1:m2_lvl2])
  Z_joint_lvl2 = t(joint_score_lvl2) %*% joint_score_lvl2 / (N*R-1)
  Z_joint_lvl2.eigen = eigen(Z_joint_lvl2)
  evalue_lvl2_multi_est = Z_joint_lvl2.eigen$values
  eigen_lvl2_var1_multi_est = mod_var1$efunctions$level2[,1:m1_lvl2] %*%
    Z_joint_lvl2.eigen$vectors[1:m1_lvl2, ]
  eigen_lvl2_var2_multi_est = mod_var2$efunctions$level2[,1:m2_lvl2] %*%
    Z_joint_lvl2.eigen$vectors[c(m1_lvl2 + 1:m2_lvl2), ]
  score_lvl2_multi_est = mod_var1$scores$level2[,1:m1_lvl2] %*%
    Z_joint_lvl2.eigen$vectors[1:m1_lvl2, ] +
    mod_var2$scores$level2[,1:m2_lvl2] %*%
    Z_joint_lvl2.eigen$vectors[c(m1_lvl2 + 1:m2_lvl2), ]
  # make prediction using univariate MFPCA results
  m_uni_pred_lvl1_var1 = min(which(mufpca_pve*cumsum(mod_var1$evalues$level1)/sum(mod_var1$evalues$level1) >= pred_pve))
  m_uni_pred_lvl1_var2 = min(which(mufpca_pve*cumsum(mod_var2$evalues$level1)/sum(mod_var2$evalues$level1) >= pred_pve))
  m_uni_pred_lvl2_var1 = min(which(mufpca_pve*cumsum(mod_var1$evalues$level2)/sum(mod_var1$evalues$level2) >= pred_pve))
  m_uni_pred_lvl2_var2 = min(which(mufpca_pve*cumsum(mod_var2$evalues$level2)/sum(mod_var2$evalues$level2) >= pred_pve))
  # make prediction using multivariate MFPCA results
  m_multi_pred_lvl1 = min(which(mufpca_pve*cumsum(evalue_lvl1_multi_est)/sum(evalue_lvl1_multi_est) >= pred_pve))
  m_multi_pred_lvl2 = min(which(mufpca_pve*cumsum(evalue_lvl2_multi_est)/sum(evalue_lvl2_multi_est) >= pred_pve))
  # run only one loop to save time
  for (i in 1:N){
    for (r in 1:R){
      index = R*(i-1) + r
      # univariate MFPCA
      z_uni_pred_var1[index,] = mod_var1$mu + 
        mod_var1$scores$level1[i,c(1:m_uni_pred_lvl1_var1)] %*%
        t(mod_var1$efunctions$level1[,c(1:m_uni_pred_lvl1_var1)]) +
        mod_var1$scores$level2[index,c(1:m_uni_pred_lvl2_var1)] %*% 
        t(mod_var1$efunctions$level2[,c(1:m_uni_pred_lvl2_var1)])
      z_uni_pred_var2[index,] = mod_var2$mu + 
        mod_var2$scores$level1[i,c(1:m_uni_pred_lvl1_var2)] %*% 
        t(mod_var2$efunctions$level1[,c(1:m_uni_pred_lvl1_var2)]) +
        mod_var2$scores$level2[index,c(1:m_uni_pred_lvl2_var2)] %*% 
        t(mod_var2$efunctions$level2[,c(1:m_uni_pred_lvl2_var2)])
      # multivariate MFPCA
      z_multi_pred_var1[index,] = mod_var1$mu +
        score_lvl1_multi_est[i, c(1:m_multi_pred_lvl1)] %*% 
        t(eigen_lvl1_var1_multi_est[,c(1:m_multi_pred_lvl1)]) +
        score_lvl2_multi_est[index, c(1:m_multi_pred_lvl2)] %*%
        t(eigen_lvl2_var1_multi_est[,c(1:m_multi_pred_lvl2)])
      z_multi_pred_var2[index,] = mod_var2$mu +
        score_lvl1_multi_est[i, c(1:m_multi_pred_lvl1)] %*% 
        t(eigen_lvl1_var2_multi_est[,c(1:m_multi_pred_lvl1)]) +
        score_lvl2_multi_est[index, c(1:m_multi_pred_lvl2)] %*%
        t(eigen_lvl2_var2_multi_est[,c(1:m_multi_pred_lvl2)])
    }
  }
  
  z_data = list(z_var1 = z_var1,
                z_var2 = z_var2)
  uni = list(eigen_lvl1_var1_uni_est = mod_var1$efunctions$level1[,1:m_uni_pred_lvl1_var1]/sqrt(2),
             eigen_lvl1_var2_uni_est = mod_var2$efunctions$level1[,1:m_uni_pred_lvl1_var2]/sqrt(2),
             eigen_lvl2_var1_uni_est = mod_var1$efunctions$level2[,1:m_uni_pred_lvl2_var1]/sqrt(2),
             eigen_lvl2_var2_uni_est = mod_var2$efunctions$level2[,1:m_uni_pred_lvl2_var2]/sqrt(2),
             evalue_lvl1_var1_uni_est = mod_var1$evalues$level1[1:m_uni_pred_lvl1_var1]*2,
             evalue_lvl1_var2_uni_est = mod_var2$evalues$level1[1:m_uni_pred_lvl1_var2]*2,
             evalue_lvl2_var1_uni_est = mod_var1$evalues$level2[1:m_uni_pred_lvl2_var1]*2,
             evalue_lvl2_var2_uni_est = mod_var2$evalues$level2[1:m_uni_pred_lvl2_var2]*2,
             fve_lvl1_var1_uni_est = mufpca_pve*mod_var1$evalues$level1[1:m_uni_pred_lvl1_var1]/
               sum(mod_var1$evalues$level1),
             fve_lvl1_var2_uni_est = mufpca_pve*mod_var2$evalues$level1[1:m_uni_pred_lvl1_var2]/
               sum(mod_var2$evalues$level1),
             fve_lvl2_var1_uni_est = mufpca_pve*mod_var1$evalues$level2[1:m_uni_pred_lvl2_var1]/
               sum(mod_var1$evalues$level2),
             fve_lvl2_var2_uni_est = mufpca_pve*mod_var2$evalues$level2[1:m_uni_pred_lvl2_var2]/
               sum(mod_var2$evalues$level2),
             z_mu_var1_uni_est = mod_var1$mu,
             z_mu_var2_uni_est = mod_var2$mu,
             # score_lvl1_var1_uni_est = mod_var1$scores$level1,
             # score_lvl1_var2_uni_est = mod_var2$scores$level1,
             # score_lvl2_var1_uni_est = mod_var1$scores$level2,
             # score_lvl3_var2_uni_est = mod_var2$scores$level2,
             z_uni_pred_var1 = z_uni_pred_var1,
             z_uni_pred_var2 = z_uni_pred_var2,
             sigma2_e_var1 = mod_var1$sigma2,
             sigma2_e_var2 = mod_var2$sigma2)
  multi = list(eigen_lvl1_var1_multi_est = eigen_lvl1_var1_multi_est[,1:m_multi_pred_lvl1],
               eigen_lvl1_var2_multi_est = eigen_lvl1_var2_multi_est[,1:m_multi_pred_lvl1],
               eigen_lvl2_var1_multi_est = eigen_lvl2_var1_multi_est[,1:m_multi_pred_lvl2],
               eigen_lvl2_var2_multi_est = eigen_lvl2_var2_multi_est[,1:m_multi_pred_lvl2],
               evalue_lvl1_multi_est = evalue_lvl1_multi_est[1:m_multi_pred_lvl1],
               evalue_lvl2_multi_est = evalue_lvl2_multi_est[1:m_multi_pred_lvl2],
               fve_lvl1_multi_est = mufpca_pve*evalue_lvl1_multi_est[1:m_multi_pred_lvl1]/
                 sum(evalue_lvl1_multi_est),
               fve_lvl2_multi_est = mufpca_pve*evalue_lvl2_multi_est[1:m_multi_pred_lvl2]/
                 sum(evalue_lvl2_multi_est),
               # score_lvl1_multi_est = score_lvl1_multi_est,
               # score_lvl2_multi_est = score_lvl2_multi_est,
               z_multi_pred_var1 = z_multi_pred_var1,
               z_multi_pred_var2 = z_multi_pred_var2)
  result = list(z_data = z_data,
                uni = uni,
                multi = multi)
  
  return(result)
}

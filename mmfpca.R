require(Matrix)
require(MASS)
require(splines)

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




mfpca.face_center = function (Y, id, visit = NULL, twoway = TRUE, weight = "obs", 
                              argvals = NULL, pve = 0.99, npc = NULL, p = 3, m = 2, knots = 35, 
                              silent = TRUE, mu) 
{
  pspline.setting.mfpca <- function(x, knots = 35, p = 3, 
                                    m = 2, weight = NULL, type = "full", knots.option = "equally-spaced") {
    K = length(knots) - 2 * p - 1
    B = spline.des(knots = knots, x = x, ord = p + 1, outer.ok = TRUE, 
                   sparse = TRUE)$design
    bs = "ps"
    if (knots.option == "quantile") {
      bs = "bs"
    }
    s.object = s(x = x, bs = bs, k = K + p, m = c(p - 1, 
                                                  2), sp = NULL)
    object = smooth.construct(s.object, data = data.frame(x = x), 
                              knots = list(x = knots))
    P = object$S[[1]]
    if (knots.option == "quantile") 
      P = P/max(abs(P)) * 10
    if (is.null(weight)) 
      weight <- rep(1, length(x))
    if (type == "full") {
      Sig = crossprod(matrix.multiply.mfpca(B, weight, 
                                            option = 2), B)
      eSig = eigen(Sig)
      V = eSig$vectors
      E = eSig$values
      if (min(E) <= 1e-07) {
        E <- E + 1e-06
      }
      Sigi_sqrt = matrix.multiply.mfpca(V, 1/sqrt(E)) %*% 
        t(V)
      tUPU = Sigi_sqrt %*% (P %*% Sigi_sqrt)
      Esig = eigen(tUPU, symmetric = TRUE)
      U = Esig$vectors
      s = Esig$values
      s[(K + p - m + 1):(K + p)] = 0
      A = B %*% (Sigi_sqrt %*% U)
    }
    if (type == "simple") {
      A = NULL
      s = NULL
      Sigi_sqrt = NULL
      U = NULL
    }
    List = list(A = A, B = B, s = s, Sigi.sqrt = Sigi_sqrt, 
                U = U, P = P)
    return(List)
  }
  quadWeights.mfpca <- function(argvals, method = "trapezoidal") {
    ret <- switch(method, trapezoidal = {
      D <- length(argvals)
      1/2 * c(argvals[2] - argvals[1], argvals[3:D] - 
                argvals[1:(D - 2)], argvals[D] - argvals[D - 
                                                           1])
    }, midpoint = c(0, diff(argvals)), stop("function quadWeights: choose either trapezoidal or midpoint quadrature rule"))
    return(ret)
  }
  if (silent == FALSE) 
    print("Organize the input")
  stopifnot((!is.null(Y) & !is.null(id)))
  stopifnot(is.matrix(Y))
  if (!is.null(visit)) {
    visit <- as.factor(visit)
  }
  else {
    visit <- as.factor(ave(id, id, FUN = seq_along))
  }
  id <- as.factor(id)
  df <- data.frame(id = id, visit = visit, Y = I(Y))
  rm(id, visit, Y)
  J <- length(levels(df$visit))
  L <- ncol(df$Y)
  nVisits <- data.frame(table(df$id))
  colnames(nVisits) = c("id", "numVisits")
  ID = sort(unique(df$id))
  I <- length(ID)
  if (is.null(argvals)) 
    argvals <- seq(0, 1, length.out = L)
  if (silent == FALSE) 
    print("Estimate population and visit-specific mean functions")
  # meanY <- colMeans(df$Y, na.rm = TRUE)
  # fit_mu <- gam(meanY ~ s(argvals))
  # mu <- as.vector(predict(fit_mu, newdata = data.frame(argvals = argvals)))
  # rm(meanY, fit_mu)
  mueta = matrix(0, L, J)
  eta = matrix(0, L, J)
  colnames(mueta) <- colnames(eta) <- levels(df$visit)
  Ytilde <- matrix(NA, nrow = nrow(df$Y), ncol = ncol(df$Y))
  if (twoway == TRUE) {
    for (j in 1:J) {
      ind_j <- which(df$visit == levels(df$visit)[j])
      if (length(ind_j) > 1) {
        meanYj <- colMeans(df$Y[ind_j, ], na.rm = TRUE)
      }
      else {
        meanYj <- df$Y[ind_j, ]
      }
      fit_mueta <- gam(meanYj ~ s(argvals))
      mueta[, j] <- predict(fit_mueta, newdata = data.frame(argvals = argvals))
      eta[, j] <- mueta[, j] - mu
      Ytilde[ind_j, ] <- df$Y[ind_j, ] - matrix(mueta[, 
                                                      j], nrow = length(ind_j), ncol = L, byrow = TRUE)
    }
    rm(meanYj, fit_mueta, ind_j, j)
  }
  else {
    Ytilde <- df$Y - matrix(mu, nrow = nrow(df$Y), ncol = L, 
                            byrow = TRUE)
  }
  df$Ytilde <- I(Ytilde)
  rm(Ytilde)
  if (silent == FALSE) 
    print("Prepare ingredients for FACE")
  if (length(knots) == 1) {
    if (knots + p >= L) 
      cat("Too many knots!\n")
    stopifnot(knots + p < L)
    K.p <- knots
    knots <- seq(-p, K.p + p, length = K.p + 1 + 2 * p)/K.p
    knots <- knots * (max(argvals) - min(argvals)) + min(argvals)
  }
  if (length(knots) > 1) 
    K.p <- length(knots) - 2 * p - 1
  if (K.p >= L) 
    cat("Too many knots!\n")
  stopifnot(K.p < L)
  c.p <- K.p + p
  List <- pspline.setting.mfpca(argvals, knots, p, m)
  B <- List$B
  Sigi.sqrt <- List$Sigi.sqrt
  s <- List$s
  U <- List$U
  A0 <- Sigi.sqrt %*% U
  G <- crossprod(B)/nrow(B)
  eig_G <- eigen(G, symmetric = T)
  G_half <- eig_G$vectors %*% diag(sqrt(eig_G$values)) %*% 
    t(eig_G$vectors)
  G_invhalf <- eig_G$vectors %*% diag(1/sqrt(eig_G$values)) %*% 
    t(eig_G$vectors)
  Bnew <- as.matrix(B %*% G_invhalf)
  Anew <- G_half %*% A0
  rm(List, Sigi.sqrt, U, G, eig_G, G_half)
  if (silent == FALSE) 
    print("Estimate the total covariance (Kt)")
  Ji <- as.numeric(table(df$id))
  diagD <- rep(Ji, Ji)
  smooth.Gt = face.Cov.mfpca(Y = unclass(df$Ytilde), argvals, 
                             A0, B, Anew, Bnew, G_invhalf, s)
  if (sum(is.na(df$Ytilde)) > 0) {
    df$Ytilde[which(is.na(df$Ytilde))] <- smooth.Gt$Yhat[which(is.na(df$Ytilde))]
  }
  if (weight == "subj") {
    YH <- unclass(df$Ytilde) * sqrt(nrow(df$Ytilde)/(I * 
                                                       diagD))
    smooth.Gt <- face.Cov.mfpca(Y = YH, argvals, A0, B, 
                                Anew, Bnew, G_invhalf, s)
    rm(YH)
  }
  diag_Gt <- colMeans(df$Ytilde^2)
  if (silent == FALSE) 
    print("Estimate principal components of the within covariance (Kw)")
  inx_row_ls <- split(1:nrow(df$Ytilde), f = factor(df$id, 
                                                    levels = unique(df$id)))
  Ysubm <- t(vapply(inx_row_ls, function(x) colSums(df$Ytilde[x, 
                                                              , drop = FALSE], na.rm = TRUE), numeric(L)))
  if (weight == "obs") {
    weights <- sqrt(nrow(df$Ytilde)/(sum(diagD) - nrow(df$Ytilde)))
    YR <- do.call("rbind", lapply(1:I, function(x) {
      weights * sqrt(Ji[x]) * t(t(df$Ytilde[inx_row_ls[[x]], 
                                            , drop = FALSE]) - Ysubm[x, ]/Ji[x])
    }))
  }
  if (weight == "subj") {
    weights <- sqrt(nrow(df$Ytilde)/sum(Ji > 1))
    YR <- do.call("rbind", lapply(1:I, function(x) {
      if (Ji[x] > 1) 
        return((weights/sqrt(Ji[x] - 1)) * t(t(df$Ytilde[inx_row_ls[[x]], 
                                                         , drop = FALSE]) - Ysubm[x, ]/Ji[x]))
    }))
  }
  smooth.Gw <- face.Cov.mfpca(Y = YR, argvals, A0, B, Anew, 
                              Bnew, G_invhalf, s)
  sigma.Gw <- smooth.Gw$evalues
  per <- cumsum(sigma.Gw)/sum(sigma.Gw)
  N.Gw <- ifelse(is.null(npc), min(which(per > pve)), min(npc, 
                                                          length(sigma.Gw)))
  smooth.Gw$efunctions <- smooth.Gw$efunctions[, 1:N.Gw]
  smooth.Gw$evalues <- smooth.Gw$evalues[1:N.Gw]
  rm(Ji, diagD, inx_row_ls, weights, weight, Ysubm, YR, B, 
     Anew, G_invhalf, s, per, N.Gw, sigma.Gw)
  if (silent == FALSE) 
    print("Estimate principal components of the between covariance (Kb)")
  temp = smooth.Gt$decom - smooth.Gw$decom
  Eigen <- eigen(temp, symmetric = TRUE)
  Sigma <- Eigen$values
  d <- Sigma[1:c.p]
  d <- d[d > 0]
  per <- cumsum(d)/sum(d)
  N.Gb <- ifelse(is.null(npc), min(which(per > pve)), min(npc, 
                                                          length(d)))
  smooth.Gb <- list(evalues = Sigma[1:N.Gb], efunctions = Bnew %*% 
                      Eigen$vectors[, 1:N.Gb])
  rm(smooth.Gt, temp, Eigen, Sigma, d, per, N.Gb, Bnew)
  if (silent == FALSE) 
    print("Estimate eigenvalues and eigenfunctions at two levels")
  efunctions <- list(level1 = as.matrix(smooth.Gb$efunctions), 
                     level2 = as.matrix(smooth.Gw$efunctions))
  evalues <- list(level1 = smooth.Gb$evalues, level2 = smooth.Gw$evalues)
  npc <- list(level1 = length(evalues[[1]]), level2 = length(evalues[[2]]))
  names(efunctions) <- names(evalues) <- names(npc) <- c("level1", 
                                                         "level2")
  rm(smooth.Gb, smooth.Gw)
  if (silent == FALSE) 
    print("Estimate the measurement error variance (sigma^2)")
  cov.hat <- lapply(c("level1", "level2"), function(x) colSums(t(efunctions[[x]]^2) * 
                                                                 evalues[[x]]))
  T.len <- argvals[L] - argvals[1]
  T1.min <- min(which(argvals >= argvals[1] + 0.25 * T.len))
  T1.max <- max(which(argvals <= argvals[L] - 0.25 * T.len))
  DIAG <- (diag_Gt - cov.hat[[1]] - cov.hat[[2]])[T1.min:T1.max]
  w2 <- quadWeights.mfpca(argvals[T1.min:T1.max], method = "trapezoidal")
  sigma2 <- max(weighted.mean(DIAG, w = w2, na.rm = TRUE), 
                0)
  rm(cov.hat, T.len, T1.min, T1.max, DIAG, w2)
  if (silent == FALSE) 
    print("Estimate principal component scores")
  Xhat <- Xhat.subject <- matrix(0, nrow(df$Y), L)
  phi1 <- efunctions[[1]]
  phi2 <- efunctions[[2]]
  score1 <- matrix(0, I, npc[[1]])
  score2 <- matrix(0, nrow(df$Y), npc[[2]])
  unVisits <- unique(nVisits$numVisits)
  if (length(unVisits) < I) {
    for (j in 1:length(unVisits)) {
      Jm <- unVisits[j]
      if (sigma2 < 1e-04) {
        A <- Jm * (t(phi1) %*% phi1)
        B <- matrix(rep(t(phi1) %*% phi2, Jm), nrow = npc[[1]])
        temp <- ginv(t(phi2) %*% phi2)
      }
      else {
        if (length(evalues[[1]]) == 1) {
          A <- Jm * (t(phi1) %*% phi1)/sigma2 + 1/evalues[[1]]
        }
        else {
          A <- Jm * (t(phi1) %*% phi1)/sigma2 + diag(1/evalues[[1]])
        }
        B = matrix(rep(t(phi1) %*% phi2/sigma2, Jm), 
                   nrow = npc[[1]])
        if (length(evalues[[2]]) == 1) {
          temp = ginv(t(phi2) %*% phi2/sigma2 + 1/evalues[[2]])
        }
        else {
          temp = ginv(t(phi2) %*% phi2/sigma2 + diag(1/evalues[[2]]))
        }
      }
      C <- t(B)
      invD <- kronecker(diag(1, Jm), temp)
      MatE <- ginv(A - B %*% invD %*% C)
      MatF <- -invD %*% C %*% MatE
      MatG <- -MatE %*% B %*% invD
      MatH <- invD - invD %*% C %*% MatG
      Mat1 <- cbind(MatE, MatG)
      Mat2 <- cbind(MatF, MatH)
      ind.Jm <- nVisits$id[which(nVisits$numVisits == 
                                   Jm)]
      YJm <- matrix(df$Ytilde[which(df$id %in% ind.Jm), 
      ], ncol = L)
      int1 <- rowsum(df$Ytilde[which(df$id %in% ind.Jm), 
      ] %*% phi1, rep(1:length(ind.Jm), each = Jm))
      int2 <- t(matrix(t(df$Ytilde[which(df$id %in% ind.Jm), 
      ] %*% phi2), nrow = npc[[2]] * Jm))
      int <- cbind(int1, int2)
      if (sigma2 >= 1e-04) {
        int <- int/sigma2
      }
      score1[which(nVisits$id %in% ind.Jm), ] <- int %*% 
        t(Mat1)
      score2[which(df$id %in% ind.Jm), ] <- t(matrix(Mat2 %*% 
                                                       t(int), nrow = npc[[2]]))
      temp <- score1[which(nVisits$id %in% ind.Jm), ] %*% 
        t(phi1)
      Xhat.subject[which(df$id %in% ind.Jm), ] <- temp[rep(1:length(ind.Jm), 
                                                           each = Jm), ]
      Xhat[which(df$id %in% ind.Jm), ] <- Xhat.subject[which(df$id %in% 
                                                               ind.Jm), ] + score2[which(df$id %in% ind.Jm), 
                                                               ] %*% t(phi2)
    }
    for (g in 1:length(levels(df$visit))) {
      ind.visit <- which(df$visit == levels(df$visit)[g])
      Xhat.subject[ind.visit, ] <- t(t(Xhat.subject[ind.visit, 
      ]) + mu + eta[, levels(df$visit)[g]])
      Xhat[ind.visit, ] <- t(t(Xhat[ind.visit, ]) + mu + 
                               eta[, levels(df$visit)[g]])
    }
    rm(YJm, g, ind.visit, ind.Jm)
  }
  else {
    for (m in 1:I) {
      Jm <- nVisits[m, 2]
      if (sigma2 < 1e-04) {
        A <- Jm * (t(phi1) %*% phi1)
        B <- matrix(rep(t(phi1) %*% phi2, Jm), nrow = npc[[1]])
        temp <- ginv(t(phi2) %*% phi2)
      }
      else {
        if (length(evalues[[1]]) == 1) {
          A <- Jm * (t(phi1) %*% phi1)/sigma2 + 1/evalues[[1]]
        }
        else {
          A <- Jm * (t(phi1) %*% phi1)/sigma2 + diag(1/evalues[[1]])
        }
        B = matrix(rep(t(phi1) %*% phi2/sigma2, Jm), 
                   nrow = npc[[1]])
        if (length(evalues[[2]]) == 1) {
          temp = ginv(t(phi2) %*% phi2/sigma2 + 1/evalues[[2]])
        }
        else {
          temp = ginv(t(phi2) %*% phi2/sigma2 + diag(1/evalues[[2]]))
        }
      }
      C <- t(B)
      invD <- kronecker(diag(1, Jm), temp)
      MatE <- ginv(A - B %*% invD %*% C)
      MatF <- -invD %*% C %*% MatE
      MatG <- -MatE %*% B %*% invD
      MatH <- invD - invD %*% C %*% MatG
      Mat1 <- cbind(MatE, MatG)
      Mat2 <- cbind(MatF, MatH)
      int1 <- colSums(matrix(df$Ytilde[df$id == ID[m], 
      ], ncol = L) %*% phi1)
      int2 <- matrix(df$Ytilde[df$id == ID[m], ], ncol = L) %*% 
        phi2
      if (sigma2 < 1e-04) {
        int <- c(int1, as.vector(t(int2)))
      }
      else {
        int <- c(int1, as.vector(t(int2)))/sigma2
      }
      score1[m, ] <- Mat1 %*% int
      score2[which(df$id == ID[m]), ] <- matrix(Mat2 %*% 
                                                  int, ncol = npc[[2]], byrow = TRUE)
      for (j in which(df$id == ID[m])) {
        Xhat.subject[j, ] <- as.matrix(mu) + eta[, df$visit[j]] + 
          as.vector(phi1 %*% score1[m, ])
        Xhat[j, ] <- Xhat.subject[j, ] + as.vector(phi2 %*% 
                                                     score2[j, ])
      }
    }
  }
  scores <- list(level1 = score1, level2 = score2)
  rm(A, B, C, int, int1, int2, invD, Mat1, Mat2, MatE, MatF, 
     MatG, MatH, temp, j, Jm, unVisits, phi1, phi2, score1, 
     score2)
  if (silent == FALSE) 
    print("Organize the results")
  res <- list(Xhat = Xhat, Xhat.subject = Xhat.subject, Y = df$Y, 
              mu = mu, eta = eta, scores = scores, efunctions = efunctions, 
              evalues = evalues, npc = npc, sigma2 = sigma2)
  rm(df, efunctions, eta, evalues, mueta, nVisits, npc, scores, 
     Xhat, Xhat.subject, argvals, diag_Gt, I, ID, J, mu, 
     pve, L, sigma2)
  return(res)
}




## The function implements the face algorithm 
face.Cov.mfpca <- function(Y, argvals, A0, B, Anew, Bnew, G_invhalf, s, Cov=FALSE, pve=0.99, npc=NULL, lambda=NULL, alpha=0.7, 
                           search.grid=TRUE, search.length=100, lower=-20, upper=20){
  
  ######## precalculation for missing data ########
  imputation <- FALSE
  Niter.miss <- 1
  L <- ncol(Y)
  n <- nrow(Y)
  
  Index.miss <- is.na(Y)
  if(sum(Index.miss)>0){
    num.miss <- rowSums(is.na(Y))
    for(i in 1:n){
      if(num.miss[i]>0){
        y <- Y[i,]
        seq <- (1:L)[!is.na(y)]
        seq2 <-(1:L)[is.na(y)]
        t1 <- argvals[seq]
        t2 <- argvals[seq2]
        fit <- smooth.spline(t1,y[seq])
        temp <- predict(fit,t2,all.knots=TRUE)$y
        if(max(t2)>max(t1)) temp[t2>max(t1)] <- mean(y[seq])
        if(min(t2)<min(t1)) temp[t2<min(t1)] <- mean(y[seq])
        Y[i,seq2] <- temp
      }
    }
    imputation <- TRUE
    Niter.miss <- 100
  }
  convergence.vector <- rep(0,Niter.miss)
  iter.miss <- 1
  lambda.input <- lambda
  totalmiss <- mean(Index.miss)
  
  
  while(iter.miss <= Niter.miss&&convergence.vector[iter.miss]==0) {
    ###################################################
    ######## Transform the Data           #############
    ###################################################
    Ytilde <- t(as.matrix(Y%*%B) %*% A0)
    C_diag <- rowSums(Ytilde^2)
    
    ###################################################
    ########  Select Smoothing Parameters #############
    ###################################################
    Y_square <- sum(Y^2)
    Ytilde_square <- sum(Ytilde^2)
    face_gcv <- function(x) {
      lambda <- exp(x)
      lambda_s <- (lambda*s)^2/(1 + lambda*s)^2
      gcv <- sum(C_diag*lambda_s) - Ytilde_square + Y_square
      trace <- sum(1/(1+lambda*s))
      gcv <- gcv/(1-alpha*trace/L/(1-totalmiss))^2
      return(gcv)
    }
    
    
    if(is.null(lambda.input) && iter.miss<=2) {
      if(!search.grid){
        fit <- optim(0,face_gcv,lower=lower,upper=upper)
        if(fit$convergence>0) {
          expression <- paste("Smoothing failed! The code is:",fit$convergence)
          print(expression)
        }
        lambda <- exp(fit$par)
      } else {
        Lambda <- seq(lower,upper,length=search.length)
        Length <- length(Lambda)
        Gcv <- rep(0,Length)
        for(i in 1:Length)
          Gcv[i] <- face_gcv(Lambda[i])
        i0 <- which.min(Gcv)
        lambda <- exp(Lambda[i0])
      }
    }
    YS <- matrix.multiply.mfpca(Ytilde,1/(1+lambda*s),2)
    
    ###################################################
    ####  Eigendecomposition of Smoothed Data #########
    ###################################################
    temp0 <- YS%*%t(YS)/n
    temp <- as.matrix(Anew%*%as.matrix(temp0%*%t(Anew)))
    Eigen <- eigen(temp,symmetric=TRUE)
    A = Eigen$vectors
    Phi = Bnew %*% A
    Sigma = Eigen$values
    
    if(iter.miss>1&&iter.miss< Niter.miss) {
      diff <- norm(YS-YS.temp,"F")/norm(YS,"F")
      if(diff <= 0.02)
        convergence.vector[iter.miss+1] <- 1
    }
    
    YS.temp <- YS
    iter.miss <- iter.miss + 1
    N <- min(n, ncol(B))
    d <- Sigma[1:N]
    d <- d[d>0]
    per <- cumsum(d)/sum(d)
    N <- ifelse (is.null(npc), min(which(per>pve)), min(npc, length(d)))
    
    #########################################
    #######     Principal  Scores   #########
    ########   data imputation      #########
    #########################################
    if(imputation) {
      Phi.N <- Phi[,1:N, drop = FALSE]
      A.N <- G_invhalf %*% A[,1:N]
      d <- Sigma[1:N]
      sigmahat2  <-  max(mean(Y[!Index.miss]^2) -sum(Sigma),0)
      if(N>1){
        Xi <- solve(t(Phi.N)%*%Phi.N + diag(sigmahat2/d)) %*% t(as.matrix(Y%*%B) %*% A.N)
      } else{
        Xi <- solve(t(Phi.N)%*%Phi.N + sigmahat2/d) %*% t(as.matrix(Y%*%B) %*% A.N)
      }
      Yhat <- t(Phi.N %*% Xi)
      Y <- Y*(1-Index.miss) + Yhat*Index.miss
      if(sum(is.na(Y))>0) print("error")
    }
    
  } ## end of while loop
  
  Phi.N <- Phi[,1:N, drop = FALSE]
  evalues <- Sigma[1:N]
  Ktilde <- NULL
  if(Cov) {
    Ktilde <- Phi.N %*%  matrix.multiply.mfpca(t(Phi.N),evalues,2)
  }
  
  return(list(Yhat=Y, decom=temp, Ktilde=Ktilde, evalues=evalues, efunctions=Phi.N))
}

matrix.multiply.mfpca <- function(A,s,option=1){
  if(option==2)
    return(A*(s%*%t(rep(1,dim(A)[2]))))
  if(option==1)
    return(A*(rep(1,dim(A)[1])%*%t(s)))
}

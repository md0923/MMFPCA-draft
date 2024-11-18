###########################################################################################
## Description: Functions for generating multilevel multi- and uni-variate two-dimensional 
##              eigenfunctions described in the simulation section of 'Joint Modeling of 
##              Evoked and Induced Event-Related Spectral Perturbations Supplementary Materials'.
###########################################################################################
## Functions included:
## Main function:
##    1. eigenf_construct: Function constructing the multilevel multi- and uni-variate 
##                         two-dimensional eigenfunctions used in the simulation.
## Supporting functions used by main function:
##    1. pi: Function forming a two-dimensional Gaussian kernel.
###########################################################################################


eigenf_construct = function(){
  #########################################################################################
  ## Description: Function constructing the multilevel multi- and uni-variate two-dimensional 
  ##              eigenfunctions described in the simulation section of the paper.
  ## args:        none
  ## Returns:     multi: a list of multilevel multivariate two-dimensional eigenfunctions
  ##              uni: a list of multilevel univariate two-dimensional eigenfunctions
  #########################################################################################
  
  # set the grids for the x and y axis of the area of definition [0,2]x[0,1]
  t = seq(0.02,2,length=100)
  f = seq(0.02,1,length=50)
  
  # set the subject-level shape parameters
  mu_11_1 = matrix(c(0.25, 0.8), ncol=1)
  omega_11_1_inv = solve(matrix(c(0.1,0.05,0.05,0.1), ncol=2))
  mu_11_2 = matrix(c(1, 0.5), ncol=1)
  omega_11_2_inv = solve(matrix(c(0.15,0.075,0.075,0.15), ncol=2))
  mu_11_3 = matrix(c(1.75, 0.2), ncol=1)
  omega_11_3_inv = solve(matrix(c(0.1,-0.05,-0.05,0.1), ncol=2))
  
  mu_12_1 = matrix(c(0.8, 0.8), ncol=1)
  omega_12_1_inv = solve(matrix(c(0.1,-0.05,-0.05,0.12), ncol=2))
  mu_12_2 = matrix(c(1.6, 0.2), ncol=1)
  omega_12_2_inv = solve(matrix(c(0.2,-0.1,-0.1,0.24), ncol=2))
  
  mu_13_1 = matrix(c(0.5, 0.8), ncol=1)
  omega_13_1_inv = solve(matrix(c(0.15,0.04,0.04,0.08), ncol=2))
  mu_13_2 = matrix(c(1.25, 0.2), ncol=1)
  omega_13_2_inv = solve(matrix(c(0.3,0.08,0.08,0.16), ncol=2))
  
  mu_14_1 = matrix(c(0.2, 0.1), ncol=1)
  omega_14_1_inv = solve(matrix(c(0.3,0.03,0.03,0.15), ncol=2))
  mu_14_2 = matrix(c(0.8, 0.4), ncol=1)
  omega_14_2_inv = solve(matrix(c(0.2,0.02,0.02,0.1), ncol=2))
  mu_14_3 = matrix(c(1.4, 0.1), ncol=1)
  omega_14_3_inv = solve(matrix(c(0.1,0.01,0.01,0.05), ncol=2))
  
  mu_15_1 = matrix(c(0.1, 0.9), ncol=1)
  omega_15_1_inv = solve(matrix(c(0.2,0.1,0.1,0.2), ncol=2))
  mu_15_2 = matrix(c(0.6, 0.4), ncol=1)
  omega_15_2_inv = solve(matrix(c(0.15,0.075,0.075,0.15), ncol=2))
  mu_15_3 = matrix(c(1.1, 0.4), ncol=1)
  omega_15_3_inv = solve(matrix(c(0.1,0.05,0.05,0.1), ncol=2))
  mu_15_4 = matrix(c(1.6, 0.9), ncol=1)
  omega_15_4_inv = solve(matrix(c(0.05,0.025,0.025,0.05), ncol=2))
  
  # set the trial-level shape parameters
  mu_21_1 = matrix(c(0.2, 0.25), ncol=1)
  omega_21_1_inv = solve(matrix(c(0.15,0.075,0.075,0.15), ncol=2))
  mu_21_2 = matrix(c(1, 0.5), ncol=1)
  omega_21_2_inv = solve(matrix(c(0.1,0.05,0.05,0.1), ncol=2))
  mu_21_3 = matrix(c(1.8, 0.75), ncol=1)
  omega_21_3_inv = solve(matrix(c(0.15,-0.075,-0.075,0.15), ncol=2))
  
  mu_22_1 = matrix(c(0.8, 0.3), ncol=1)
  omega_22_1_inv = solve(matrix(c(0.15,-0.04,-0.04,0.08), ncol=2))
  mu_22_2 = matrix(c(1.6, 0.7), ncol=1)
  omega_22_2_inv = solve(matrix(c(0.3,-0.08,-0.08,0.16), ncol=2))
  
  mu_23_1 = matrix(c(0.4, 0.7), ncol=1)
  omega_23_1_inv = solve(matrix(c(0.1,0.05,0.05,0.12), ncol=2))
  mu_23_2 = matrix(c(1.2, 0.3), ncol=1)
  omega_23_2_inv = solve(matrix(c(0.2,0.1,0.1,0.24), ncol=2))
  
  mu_24_1 = matrix(c(0.2, 0.9), ncol=1)
  omega_24_1_inv = solve(matrix(c(0.25,-0.05,-0.05,0.5), ncol=2))
  mu_24_2 = matrix(c(0.7, 0.4), ncol=1)
  omega_24_2_inv = solve(matrix(c(0.15,-0.03,-0.03,0.3), ncol=2))
  mu_24_3 = matrix(c(1.2, 0.9), ncol=1)
  omega_24_3_inv = solve(matrix(c(0.05,-0.01,-0.01,0.1), ncol=2))
  
  mu_25_1 = matrix(c(0.75, 0.2), ncol=1)
  omega_25_1_inv = solve(matrix(c(0.3,0.03,0.03,0.15), ncol=2))
  mu_25_2 = matrix(c(1.25, 0.6), ncol=1)
  omega_25_2_inv = solve(matrix(c(0.2,0.02,0.02,0.1), ncol=2))
  mu_25_3 = matrix(c(1.75, 0.2), ncol=1)
  omega_25_3_inv = solve(matrix(c(0.1,0.01,0.01,0.05), ncol=2))
  
  # generate the subject-level basis functions
  b_11 = pi(t, f, mu_11_1, omega_11_1_inv) +
         pi(t, f, mu_11_2, omega_11_2_inv) -
         pi(t, f, mu_11_3, omega_11_3_inv)
  b_12 = pi(t, f, mu_12_1, omega_12_1_inv) -
         pi(t, f, mu_12_2, omega_12_2_inv)
  b_13 = pi(t, f, mu_13_1, omega_13_1_inv) +
         pi(t, f, mu_13_2, omega_13_2_inv)
  b_14 = pi(t, f, mu_14_1, omega_14_1_inv) +
         pi(t, f, mu_14_2, omega_14_2_inv) -
         pi(t, f, mu_14_3, omega_14_3_inv)
  b_15 = pi(t, f, mu_15_1, omega_15_1_inv) +
         pi(t, f, mu_15_2, omega_15_2_inv) -
         pi(t, f, mu_15_3, omega_15_3_inv) +
         pi(t, f, mu_15_4, omega_15_4_inv)
  
  # generate the trial-level basis functions
  b_21 = pi(t, f, mu_21_1, omega_21_1_inv) -
         pi(t, f, mu_21_2, omega_21_2_inv) -
         pi(t, f, mu_21_3, omega_21_3_inv)
  b_22 = pi(t, f, mu_22_1, omega_22_1_inv) -
         pi(t, f, mu_22_2, omega_22_2_inv)
  b_23 = pi(t, f, mu_23_1, omega_23_1_inv) +
         pi(t, f, mu_23_2, omega_23_2_inv)
  b_24 = pi(t, f, mu_24_1, omega_24_1_inv) -
         pi(t, f, mu_24_2, omega_24_2_inv) +
         pi(t, f, mu_24_3, omega_24_3_inv)
  b_25 = pi(t, f, mu_25_1, omega_25_1_inv) -
         pi(t, f, mu_25_2, omega_25_2_inv) +
         pi(t, f, mu_25_3, omega_25_3_inv)
  
  # save the generated basis functions at each level into columns of a matrix
  basis_lvl1 = matrix(c(b_11, b_12, b_13, b_14, b_15), ncol = 5)
  basis_lvl2 = matrix(c(b_21, b_22, b_23, b_24, b_25), ncol = 5)
  
  # multivariate eigenfunctions
  # gram-schmidt process
  require("pracma") # package containing gramSchmidt functions
  basis_lvl1_gs = gramSchmidt(basis_lvl1)$Q *sqrt(5000/2)
  basis_lvl2_gs = gramSchmidt(basis_lvl2)$Q *sqrt(5000/2)
  # split the area of definition into two
  eigen_lvl1_var1 = basis_lvl1_gs[1:2500,]
  eigen_lvl1_var2 = basis_lvl1_gs[2501:5000,]
  eigen_lvl2_var1 = basis_lvl2_gs[1:2500,]
  eigen_lvl2_var2 = basis_lvl2_gs[2501:5000,]
  
  # univariate eigenfunctions
  # split the area of definition into two and then apply gram-schmidt process
  eigen_lvl1_var1_single = gramSchmidt(basis_lvl1[1:2500,])$Q * sqrt(2500)
  eigen_lvl1_var2_single = gramSchmidt(basis_lvl1[2501:5000,])$Q * sqrt(2500)
  eigen_lvl2_var1_single = gramSchmidt(basis_lvl2[1:2500,])$Q * sqrt(2500)
  eigen_lvl2_var2_single = gramSchmidt(basis_lvl2[2501:5000,])$Q * sqrt(2500)
  
  # save the eigenfunctions into lists
  multi = uni = list()
  # list for the multivariate eigenfunctions
  multi$lvl1$var1 = eigen_lvl1_var1
  multi$lvl1$var2 = eigen_lvl1_var2
  multi$lvl2$var1 = eigen_lvl2_var1
  multi$lvl2$var2 = eigen_lvl2_var2
  # list for the univariate eigenfunctions
  uni$lvl1$var1 = eigen_lvl1_var1_single
  uni$lvl1$var2 = eigen_lvl1_var2_single
  uni$lvl2$var1 = eigen_lvl2_var1_single
  uni$lvl2$var2 = eigen_lvl2_var2_single
  
  result = list(multi = multi,
                uni = uni)
  return(result)
}




pi = function(
    t_arg,                         # sampling grid for T (vector of length m_t)
    f_arg,                         # sampling grid for F (vector of length m_f)
    mu,                            # mean vector of the Gaussian kernel (vector of length 2)
    omega_inv){                    # inverse of the covariance matrix of the Gaussian kernel 
                                   #  (matrix of 2x2)
  #########################################################################################
  ## Description: Function that forms a two-dimensional Gaussian kernel specified by shape
  ##              parameters (mean and covariance matrix) and evaluates at the given two- 
  ##              dimensional area of definition for TxF.
  ## args:        see above
  ## Returns:     pi_tf: the Gaussian kernel specified by the given shape parameters evaluated
  ##                     at the area of definition spanned by sampling grid for T and F 
  ##                     (matrix of m_f x m_t)
  #########################################################################################  
  
  # obtain the number of grids on each axis of the two-dimensional area of definition
  m_t = length(t_arg)
  m_f = length(f_arg)
  
  # create the result matrix that is of same dimension as the area of definition
  pi_tf = matrix(NA, nrow = m_f, ncol = m_t)
  
  # make sure the format of mean vector is compatible with matrix multiplication
  mu = matrix(c(mu), nrow = 2, ncol = 1)
  
  # evaluate the Gaussian kernel at each sampling point from the area of definition
  for(i in 1:m_f){
    for(j in 1:m_t){
      # calculate the value of Gaussian kernel at each sampling point
      tf_vec = matrix(c(t_arg[j], f_arg[i]), ncol=1)
      pi_tf[i,j] = exp(-0.5*t(tf_vec - mu) %*% omega_inv %*% (tf_vec - mu))
    }
  }
  
  return(pi_tf)
}

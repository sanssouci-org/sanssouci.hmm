
#' Title
#'
#' @param x vector of statistics
#' @param fw_bc result of 
#' @param seuil threshold for p values
#' @param A_est the estimation of the transition matrix
#' @param f0x_est estimation of the null density
#' @param f1x_est estimation of the alternative density
#' @param Pi_est estimation of the initial distribution
#' @param min_size keeping only a given number of 
#' @param min_jump 
#' @param pval a vector of pvalue
#' @param all does the set selecting all the hypothesis have to be add ?
#'
#' @return
#' @export
#'
#' @examples
Selection_delta <- function(x, fw_bc, seuil, A_est,
                            f0x_est, f1x_est, Pi_est, min_size,
                            min_jump = NULL,
                            pval = NULL, all = FALSE){
  
  m <- length(x)
  
  LIS <- enframe(fw_bc$gamma[, 1]) %>%
    arrange(value) %>%
    rowid_to_column() %>%
    mutate(k = cumsum(value) / rowid) %>%
    filter(k <= seuil) %>%
    arrange(name)
  
  sel_sc <- LIS$name
  
  
  viterbi <- viterbi_log(m, log(A_est), log(f0x_est),
                         log(f1x_est), log(Pi_est))
  
  
  sel_viter_min_size <- long_reg(viterbi, min_size)
  
  sel_viter_est <- which(viterbi == 1)
  if(is.null(pval)){
    pval <- pval_m(x, fw_bc$gamma[,1])
  }
  pvalm_tresh <- which(pval < seuil)
  pvalm_sel <- pvalm_tresh[pvalm_tresh < 200]
  
  
  Sel <- tibble( Sel = list(
    pvalm_sel,
    pvalm_tresh,
    sel_sc,
    sel_viter_est,
    sel_viter_min_size),
    Nom = c(
      "pvalm_sel",
      "pval_tresh",
      "lfdr_tresh",
      "sel_viter_est",
      "sel_viter_min_size")) %>%
    mutate(Size = map_dbl(Sel,~length(.)))
}





#' Title
#'
#' @param m numeric, the number of hypothesis
#' @param A matrix, the transition matrix of the markov chain
#' @param Pi a numeric vector, the initial distribution (if different from the stationary distribution)
#' @param SNR numeric, the Signal to noise ratio
#' @param type_sim character, the type of simulation made 
#' @param al  numeric (between 0 and 1),  the risk
#' @param seuil  numeric, the threshold for the pvalues 
#' @param h  numeric, the window in kernel density estimation
#' @param n_boot numeric,  the number of bootstrap samples
#' @param min_size numeric,  the minimum size of a block of consecutive selected hypothesis 
#' @param norm logical, whether densities are gaussian or not (otherwise it is assume to be a student)
#' @param m0  numeric, the expectency under the null
#' @param sd0 numeric,  the standard deviation under the null
#' @param df  numeric, the degres of freedom 
#' @param m0_init  numeric, if f0 is unknown the initialisation of its expectency
#' @param sd0_init numeric, if f0 is unknown the initialisation of its expectency
#' @param df_init numeric, if f0 is unknown the initialisation of its expectency
#' @param norm_init logical, if f0 is unknown wether the initialisation of its law is gaussian
#' @param max_pi0 numeric, maximum value of pi0
#' @param type_init 
#' @param num_seed to set a seed. 
#' @param f0_known logical, wether f0 is assumed to be known or not
#' @param approx logical wether to estimate the density at every point or to add a small linear approximation between a large range of point
#' @param all logical, compare all the 
#' @param size_b0 numerical, size of continuous block for type_sim = "block"
#' @param pct_b1 numerical, ratio of the size of the block une H1 / size_b0
#' @param include_H0 logical, wether to add the information on H0 or not
#' @param delta numeric, the value of delta (the part of the risk for the bootstrap and the other one)
#'
#' @return
#' @export
#'
#' @examples
#' simu_delta( m = c(100),
#' A = matrix(c(0.95, 0.05, 0.2, 0.80), 2, 2, byrow = T),
#' Pi = c(0.95, 0.05),
#' rho = c(0),
#' SNR = 2,
#' prob = c(0.5),
#' type_sim = c("HMM"),
#' n_boot = 20,
#' al = 0.2, s_dbnr = 10, b_act = 2, d = 1, seuil= 0.05,
#' min_size = 2, norm = TRUE, sd0 = 0.5, m0= 0,sd0_init = 0.5, m0_init= 0,
#' norm_init = TRUE, df= 2, num_seed= 1234, type_init="given", f0_known=TRUE,
#' approx = TRUE, delta = 0.9)
# simu_delta <- function(m, A, Pi = NULL,  rho, SNR, prob, type_sim = "HMM", al, s_dbnr,
#                        b_act, d, seuil,
#                        h =0.3,  n_boot, min_size, norm, m0, sd0, df,
#                        m0_init, sd0_init, df_init, norm_init, max_pi0= 0.99999,
#                        type_init, num_seed, f0_known, approx, all =FALSE,
#                        size_b0= 300, pct_b1 =1/3, include_H0 = FALSE, delta, n_seg, 
#                        drop1, drop2, tumorFraction, sel_function = Selection_delta) {
#   set.seed(num_seed)
#   if(is.null(Pi)){
#     Pi0 <- A[2,1] / (1 + A[2,1] - A[1,1])
#     Pi <- c(Pi0, 1- Pi0)
#   }
#   ## Simuation des donnees
#   if(type_sim =="HMM"){
#     theta <- sim_markov(m, Pi, A)
#   }
#   if (type_sim == "HMM_nonstat") {
#     theta <- sim_markov_nonstat(m, Pi)
#   }
#   
#   if(type_sim =="block"){
#     theta <- rep(rep(0:1, c(size_b0, size_b0 * pct_b1)),
#                  m / (size_b0 * (1 + pct_b1)))
#     nb0 <-  m / (1+pct_b1)
#     nb1 <- m - nb0
#     A <- matrix(c((nb0-8)/(nb0-1), (8)/(nb0-1),(7)/(nb1-1),(nb1-7)/(nb1-1)), ncol =2, byrow = TRUE)
#   }
#   if (type_sim == "realistic") {
#     realistic <-
#       sim_realistic(len = m,
#                     n_seg = n_seg,
#                     drop1,
#                     drop2,
#                     tumorFraction)
#     theta <- realistic$theta
#     x <- realistic$w
#     rm(realistic)
#     gc()
#     thet_deb <- theta[-m]
#     thet_fin <- theta[-1]
#     nb00 <- sum((thet_deb - thet_fin + 1) * (1 - thet_deb))
#     nb11 <- sum(-(thet_deb - thet_fin - 1) * (thet_deb))
#     nb0  <- sum(1 - theta)
#     nb1 <- sum(theta)
#     A <-
#       matrix(
#         c(nb00 / nb0, 1 -  nb00 / nb0, 1 - nb11 / nb1,  nb11 / nb1),
#         byrow = TRUE,
#         ncol = 2
#       )
#   }
#   if (type_sim != "realistic") {
#     x <- rep(0, m)
#     if (norm) {
#       x[theta == 0] <- rnorm(sum(theta == 0), m0, sd0)
#       x[theta == 1] <- rnorm(sum(theta == 1), SNR * sd0 + m0,
#                              sd0)
#     }
#     else {
#       x[theta == 0] <- rt(sum(theta == 0), df)
#       x[theta == 1] <- rt(sum(theta == 1), df) + SNR
#     }
#   }
#   
#   if(type_init == "locfdr"){
#     w <- locfdr(x)
#     m0_init <- w$fp0[3,1]
#     sd0_init <-  w$fp0[3,2]
#     norm_init <- TRUE
#   }
#   
#   if(norm_init){
#     pval <- 2 * (1 - pnorm(abs(x), m0_init, sd0_init))
#   }else{
#     pval <- 2 * (1 - pt(abs(x), df_init))
#   }
#   
#   
#   ## Pour HMM oracle
#   if(norm){
#     f0x <- dnorm(x, m0, sd0)
#     f1x <- dnorm(x, SNR*sd0 + m0, sd0)
#   }else{
#     f0x <- dt(x, df)
#     f1x <- dt(x-SNR, df)
#   }
#   
#   fw_bc_or <- for_back(m, A, f0x, f1x, Pi)
#   Pis_or <- lapply(2:m, function(i){
#     get_A( m,alpha = fw_bc_or$alpha, beta = fw_bc_or$beta, A, f0x, f1x, i = i)
#   } )
#   
#   
#   
#   
#  
# 
#   Final <- res_all (x, al = al, sel_function, delta = delta, h = h,
#                           norm_init = norm_init,
#                           m0_init = m0_init, sd0_init = sd0_init, 
#                           df_init = df_init,
#                           max_pi0= max_pi0,
#                           min_size = min_size,
#                           f0_known = f0_known, 
#                           approx = approx, drop_sel = FALSE)   %>%
#     mutate(
#       IC_or = map(Sel,~get_IC(sel = ., li0 = fw_bc_or$gamma[,1],
#                               Pis = Pis_or, f0x =f0x ,
#                               f1x= f1x, alpha = al)),
#       V_HMM_small_or =  map_dbl(IC_or,~.[1]),
#       V_HMM_or =  map_dbl(IC_or,~.[2]),
#       FDR_or  = map_dbl(Sel, ~sum(fw_bc_or$gamma[.,1])),
#       V_real = map_dbl(Sel , ~sum(theta[.] == 0))
#     ) %>%
#     select(-Sel,  -IC_or) 
#   return(Final)
# } 






#' Title
#' 
#' @useDynLib hmm.sanssouci
#' @importFrom Rcpp evalCpp
#' @import sansSouci
#' @importFrom  magrittr "%>%"
#' @importFrom dplyr arrange
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom tibble rowid_to_column
#' @importFrom tibble enframe
#' @importFrom tibble tibble
#' @param x numeric vector, of statistics (order using a given order)
#' @param al numeric, the risk
#' @param sel_function a function that return a tibble with selected set (see details)
#' @param delta anumeric,  risk delta to share the risk alpha between the bootstrap part and the estimated part
#' @param h numeric, the window size for the kde 
#' @param f0_known logical, wether f0 is known (if f0_known =TRUE the initialisation will be the true law under H0)
#' @param norm_init logical, wether the initialisation is normal or not (consider student)
#' @param n_boot numeric,  number of bootstrap sample
#' @param seuil numeric,  threshold for selected pvalues
#' @param m0_init numeric,  expectency under H0 (if norm_init = TRUE)
#' @param sd0_init numeric,  expectency under H0
#' @param df_init numeric, if norm_init = FALSE, student degree of freedom
#' @param max_pi0 the maximum value for the first estimation of pi0
#' @param approx wheter the kde is approximated using linear interpolation with a large range of values or calculate for every point.
#'
#' @return
#' @export
#'
#' @examples
#' m <- 2000
#' theta <- sim_markov(m, Pi = c(0.8,0.2), A = matrix(c(0.95, 0.05, 0.2, 0.80), 2, 2, byrow = T))
#' x <- rep(0, m)
#'  x[theta == 0] <- rnorm(sum(theta ==0))
#'  x[theta == 1] <- rnorm(sum(theta ==1), 2, 1)
#'  sanssouci.hmm(x, al= 0.1, sel_function = Selection_delta) 
sanssouci.hmm <- function(x, al, sel_function, delta = 0.1, h = 0.3,
                          f0_known = TRUE, 
                          norm_init = TRUE,
                          n_boot = 20,
                          seuil= 0.05,
                          m0_init = 0, sd0_init = 1, df_init = NULL,
                          max_pi0= 0.99999,
                          approx = TRUE, 
                          min_size =5, 
                          min_jump = 3, type_init = "given", drop_sel = TRUE){
  
  m <- length(x)
  ## Pour HMM est
  
  Est <- Estimation( x, h =h,
                     m0_init, sd0_init, df_init, norm_init, max_pi0= max_pi0,
                     f0_known = f0_known, f0x_est = NULL, pval = NULL,
                     plot = FALSE, size_plot = min(10000, length(x)),
                     approx = approx)
  if(norm_init){
    pval <- 2 * (1 - pnorm(abs(x), m0_init, sd0_init))
  }else{
    pval <- 2 * (1 - pt(abs(x), df_init))
  }
  
  
  if(f0_known){
    pval_don <- pval
  }else{
    pval_don = NULL
  }
  A_est <- Est$Em$A
  Pi_est <- Est$Em$Pi
  Sel <- sel_function(x, Est$Em$fw_bc_EM, seuil, A_est,
                         f0x_est =Est$Em$f0x ,
                         f1x_est= Est$Em$f1x, Pi_est, min_size,
                         pval = pval_don,
                         all = all )
  
  Pis_est <- lapply(2:m, function(i){
    get_A( m,alpha = Est$Em$fw_bc_EM$alpha, beta = Est$Em$fw_bc_EM$beta,
           A_est, f0x =Est$Em$f0x ,
           f1x= Est$Em$f1x, i = i)
  })
  
  
  ### Pour Estimer bootstrap :
  boots <- enframe( x = 1:n_boot, name = NULL, value = "id_boot") %>%
    mutate(HMM_boot = map(id_boot, ~boots_delta(A_est = A_est, Pi_est= Pi_est,
                                                x_from = x,
                                                prob1 = Est$Em$fw_bc_EM$gamma[,2], h =h,
                                                Sel = Sel, al = 0.1,
                                                seuil = seuil,
                                                min_size =  min_size,
                                                max_pi0 = max_pi0,
                                                m0_init =  m0_init, sd0_init = sd0_init,
                                                df_init  = df_init, norm_init = norm_init,
                                                type_init = type_init,
                                                approx = approx,
                                                delta = delta, 
                                                sel_function = sel_function))) %>%
    unnest(HMM_boot) %>%
    select(- Sel) %>%
    nest(Est_HMM_boot = c(id_boot, Real_boot,
                          V_HMM_est_boot_aldemi,
                          V_HMM_small_est_boot_aldemi,
                          V_HMM_oracle_boot_aldemi,
                          V_HMM_small_oracle_boot_aldemi,
                          V_HMM_est_boot_aldemi_samesel,
                          V_HMM_small_est_boot_aldemi_samesel,
                          V_HMM_oracle_boot, V_HMM_small_oracle_boot,
                          V_HMM_oracle_boot_al1_moins_delta,
                          V_HMM_small_oracle_boot_al1_moins_delta,
                          V_HMM_oracle_boot_aldelta,
                          V_HMM_small_oracle_boot_aldelta,
                          V_HMM_est_boot_aldelta,
                          V_HMM_small_est_boot_aldelta,
                          V_HMM_est_boot_al1_moins_delta,
                          V_HMM_small_est_boot_al1_moins_delta,
                          V_HMM_est_boot,
                          V_HMM_small_est_boot,
                          V_HMM_est_boot_aldelta_samesel,
                          V_HMM_small_est_boot_aldelta_samesel,
                          V_HMM_est_boot_al1_moins_delta_samesel,
                          V_HMM_small_est_boot_al1_moins_delta_samesel,
                          V_HMM_est_boot_samesel,
                          V_HMM_small_est_boot_samesel,  
                          Espe,
                          Size,
                          Size_boot,
                          Sel_from))
  
  Det_a <- det(A_est)
  
  Final <- Sel   %>%
    mutate(
      Det_A_est = Det_a,
      Borne_est = map(Sel,~get_probs(sel = ., li0 = Est$Em$fw_bc_EM$gamma[,1],
                                     Pis = Pis_est, f0x =Est$Em$f0x ,
                                     f1x= Est$Em$f1x, probs = c(al/2, 1-al/2,
                                                                al, 1-al,
                                                                al*(1-delta), 1 - al*(1-delta), 
                                                                al*(delta), 1 - al*(delta)))),
      V_simes = map_dbl(Sel,~borne ( type_borne = "Simes", sel = ., m = m,
                                     pval = pval, alpha = al)
      ),
      V_HMM_est_aldemi =  map_dbl(Borne_est,~.[2]),
      V_HMM_small_est_aldemi =  map_dbl(Borne_est,~.[1]),
      V_HMM_est_aldelta =  map_dbl(Borne_est,~.[8]),
      V_HMM_small_est_aldelta =  map_dbl(Borne_est,~.[7]),
      V_HMM_est_al1_moins_delta =  map_dbl(Borne_est,~.[6]),
      V_HMM_small_est_al1_moins_delta =  map_dbl(Borne_est,~.[5]),
      V_HMM_small_est =  map_dbl(Borne_est,~.[3]),
      V_HMM_est =  map_dbl(Borne_est,~.[4]),
      FDR_est = map_dbl(Sel, ~sum( Est$Em$fw_bc_EM$gamma[.,1]))
    ) 
  
  if(drop_sel){
    Final <- Final   %>%
      select(-Sel,  -Borne_est) %>%
      left_join(boots, by = c("Nom"))
  } else {
    Final <- Final   %>%
      select(  -Borne_est) %>%
      left_join(boots, by = c("Nom"))
    
  }

  
  return(Final)
}


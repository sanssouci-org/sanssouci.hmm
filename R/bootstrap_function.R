boots_delta <- function (A_est, Pi_est =NULL, x_from, prob1, h, Sel_from, al, seuil,
                         min_size, max_pi0, m0_init, sd0_init, df_init, norm_init,
                         type_init, approx, maxit=100,f0_known = TRUE,delta,
                         type_check = "closer_to0", 
                         sel_function = Selection_delta)
{
    Pi0 <- A_est[2,1] / (1 + A_est[2,1] - A_est[1,1])
    Pi_est <- c(Pi0, 1- Pi0)
  
  Sel_from <- Sel_from %>% rename(Sel_from = Sel)
  m <- length(x_from)
  Data_temp <- sim_hmm_from_weightkde(A_est, Pi_est, x_from,
                                      prob1, h, n)
  theta <- Data_temp$theta
  x <- Data_temp$x
  rm(Data_temp)
  gc()
  if (norm_init) {
    pval <- 2 * (1 - pnorm(abs(x), m0_init, sd0_init))
  }
  else {
    pval <- 2 * (1 - pt(abs(x), df_init))
  }
  if (type_init == "locfdr") {
    w <- locfdr(x)
    m0_init <- w$fp0[3, 1]
    sd0_init <- w$fp0[3, 2]
    norm_init <- TRUE
  }
  if (norm_init) {
    f0x <- dnorm(x, m0_init, sd0_init)
    if(f0_known){
      f0x_from <- dnorm(x_from, m0_init, sd0_init)
    }
    
  }
  else {
    f0x <- dt(x, df_init)
    if(f0_known){
      f0x_from <- dt(x_from, df_init)
    }
  }
  
  if (approx) {
    d1 <- density(x)
    
    f1x_first <- sapply(d1$x, function(xi) {
      sum(K((x_from - xi)/h) * prob1)/sum(h * prob1)
    })
    f1x <- approx(d1$x, f1x_first, x)$y
    rm(f1x_first)
    if(!f0_known){
      f0x_first <- sapply(d1$x, function(xi) {
        sum(K((x_from - xi)/h) * (1 - prob1))/sum(h * (1 -
                                                         prob1))
      })
      f0x <- approx(d1$x, f0x_first, x)$y
      rm(f0x_first)
    }
    gc()
  }
  else {
    f1x <- sapply(x, function(xi) {
      sum(K((x_from - xi)/h) * prob1)/sum(h * prob1)
    })
    if(!f0_known){
      f0x <- sapply(x, function(xi) {
        sum(K((x_from - xi)/h) * prob1)/sum(h * (1-prob1))
      })
    }
  }
  fw_bc_or_star <- for_back(m, A_est, f0x, f1x, Pi_est)
  Pis_or_star <- lapply(2:m, function(i) {
    get_A(m, alpha = fw_bc_or_star$alpha, beta = fw_bc_or_star$beta,
          A_est, f0x, f1x, i = i)
  })
  
  Est <- Estimation( x, h =h,
                     m0_init, sd0_init, df_init, norm_init, max_pi0= max_pi0,
                     f0_known = f0_known, f0x_est = NULL, pval = NULL,
                     plot = FALSE, size_plot= min(10000, length(x)),
                     approx = approx, maxit=maxit)
  
  if(f0_known){
    pval_don <- pval
  }else{
    pval_don = NULL
  }
  if(!f0_known && type_check == "closer_to0" &&
     abs(sum(x*Est$Em$f0x)) > abs(sum(x*Est$Em$f1x))){
    Est$Em$A <- Est$Em$A[2:1,2:1]
    Est$Em$Pi <-Est$Em$Pi[2:1]
    f0x_est <- Est$Em$f1x
    Est$Em$f1x <- Est$Em$f0x
    Est$Em$f0x <- f0x_est
    Est$Em$fw_bc_EM <- for_back(m, Est$Em$A ,Est$Em$f0x,
                                Est$Em$f1x, Est$Em$Pi)
    rm(f0x_est)
    gc()
  }
  A_est <- Est$Em$A
  Pi_est <- Est$Em$Pi
  browser()
  Sel <- sel_function(x, Est$Em$fw_bc_EM, seuil, A_est,
                         f0x_est =Est$Em$f0x ,
                         f1x_est= Est$Em$f1x, Pi_est, min_size,
                         pval =  pval_don,
                         all =  all ) %>%
    rename(Size_boot = Size)
  Pis_est_star <- lapply(2:m, function(i) {
    get_A(m, alpha = Est$Em$fw_bc_EM$alpha, beta = Est$Em$fw_bc_EM$beta,
          A_est, Est$Em$f0x, Est$Em$f1x, i = i)
  })
  if (approx) {
    d1_from <- density(x_from)
    f1x_first <- sapply(d1_from$x, function(xi) {
      sum(K((x - xi)/h) * Est$Em$fw_bc_EM$gamma[, 2])/sum(h * Est$Em$fw_bc_EM$gamma[, 2])
    })
    f1x_from <- approx(d1_from$x, f1x_first, x_from)$y
    if(!f0_known){
      f0x_first <- sapply(d1_from$x, function(xi) {
        sum(K((x - xi)/h) * Est$Em$fw_bc_EM$gamma[, 1])/sum(h * Est$Em$fw_bc_EM$gamma[, 1])
      })
      f0x_from <- approx(d1_from$x, f0x_first, x_from)$y 
    }
    
    rm(f1x_first)
    rm(f0x_first)
    gc()
  }
  else {
    f1x_from <- sapply(x_from, function(xi) {
      sum(K((x - xi)/h) * Est$Em$fw_bc_EM$gamma[, 2])/sum(h * Est$Em$fw_bc_EM$gamma[, 2])
    })
    if(! f0_known){
      f0x_from <- sapply(x_from, function(xi) {
        sum(K((x - xi)/h) * Est$Em$fw_bc_EM$gamma[, 1])/sum(h * Est$Em$fw_bc_EM$gamma[, 1])
      })}
  }
  fw_bc_from <- for_back(m, A_est, f0x_from, f1x_from,
                         Pi_est)
  Pis_from <- lapply(2:m, function(i) {
    get_A(m, alpha = fw_bc_from$alpha, beta = fw_bc_from$beta,
          A_est, f0x_from, f1x_from, i = i)
  })
  prob <- c(al/2, 1-al/2,
            al*(1-delta), 1 - al*(1-delta),
            al, 1-al,
            al*(delta), 1 - al*(delta))
  
  Sel_boot <- Sel %>% full_join(Sel_from, by = "Nom") %>%
    mutate(Borne_from = map(Sel_from,
                            ~get_probs(sel = ., li0 = fw_bc_from$gamma[, 1],
                                       Pis = Pis_from, f0x = f0x_from, f1x = f1x_from, prob)),
           Borne_oracle_boot = map(Sel, ~get_probs(sel = .,
                                                   li0 = fw_bc_or_star$gamma[, 1], Pis = Pis_or_star,
                                                   f0x = f0x, f1x = f1x, prob)),
           Borne_est_boot = map(Sel,   ~get_probs(sel = ., li0 = Est$Em$fw_bc_EM$gamma[, 1],
                                                  Pis = Pis_est_star, f0x = Est$Em$f0x,
                                                  f1x = Est$Em$f1x, prob)),
           V_HMM_oracle_boot = map_dbl(Borne_oracle_boot, ~.[6]),
           V_HMM_small_oracle_boot = map_dbl(Borne_oracle_boot, ~.[5]),
           V_HMM_oracle_boot_aldemi = map_dbl(Borne_oracle_boot, ~.[2]),
           V_HMM_small_oracle_boot_aldemi = map_dbl(Borne_oracle_boot, ~.[1]),
           V_HMM_oracle_boot_al1_moins_delta = map_dbl(Borne_oracle_boot, ~.[4]),
           V_HMM_small_oracle_boot_al1_moins_delta = map_dbl(Borne_oracle_boot, ~.[3]),
           V_HMM_oracle_boot_aldelta = map_dbl(Borne_oracle_boot, ~.[8]),
           V_HMM_small_oracle_boot_aldelta = map_dbl(Borne_oracle_boot, ~.[7]),
           V_HMM_est_boot_aldemi = map_dbl(Borne_est_boot, ~.[2]),
           V_HMM_small_est_boot_aldemi = map_dbl(Borne_est_boot, ~.[1]),
           V_HMM_est_boot_aldelta = map_dbl(Borne_est_boot, ~.[8]),
           V_HMM_small_est_boot_aldelta = map_dbl(Borne_est_boot, ~.[7]),
           V_HMM_est_boot_al1_moins_delta = map_dbl(Borne_est_boot, ~.[4]),
           V_HMM_small_est_boot_al1_moins_delta = map_dbl(Borne_est_boot, ~.[3]),
           V_HMM_est_boot = map_dbl(Borne_est_boot, ~.[6]),
           V_HMM_small_est_boot = map_dbl(Borne_est_boot, ~.[5]),
           V_HMM_est_boot_aldemi_samesel = map_dbl(Borne_from, ~.[2]),
           V_HMM_small_est_boot_aldemi_samesel = map_dbl(Borne_from, ~.[1]),
           V_HMM_est_boot_aldelta_samesel = map_dbl(Borne_from, ~.[4]),
           V_HMM_small_est_boot_aldelta_samesel = map_dbl(Borne_from, ~.[3]),
           
           V_HMM_est_boot_al1_moins_delta_samesel = map_dbl(Borne_from, ~.[8]),
           V_HMM_small_est_boot_al1_moins_delta_samesel = map_dbl(Borne_from, ~.[7]),
           V_HMM_est_boot_samesel = map_dbl(Borne_from, ~.[6]),
           V_HMM_small_est_boot_samesel = map_dbl(Borne_from, ~.[5]),
           Espe = map_dbl(Sel, ~sum(fw_bc_or_star$gamma[.,1])),
           Real_boot = map_dbl(Sel, ~sum(theta[.] == 0)),
           Est = list(Est)) %>%
    select(-Borne_oracle_boot, -Borne_est_boot, -Borne_from)
  return(Sel_boot)
}

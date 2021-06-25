
quant_max <- function(al, base){
  a<-  base %>%
    mutate(eps =  map_dbl(eps, ~min(max(.x,-1), 1))) %>%
    pull(eps) %>% sort()
  if(al >= 0.5){
    max(a[ceiling(length(a)*al)],0)
  }else{
    min(a[floor(length(a)*al)+1],0)
  }
  
}
mise_en_forme  <-  function(Result) {
  alpha <- Result$al %>% unique()
  delt <- Result$delta %>% unique()
  Res_all <-   Result %>%
    filter(Size > 0) %>%
    mutate(
      q_b3_high = map_dbl(Est_HMM_boot, function(x) {
        base <- mutate(x,
                       eps =   (Real_boot - V_HMM_est_boot) / Size_boot)
        quant_max(1 - alpha, base)
      }),# in a perfect word q_b3_high will be 0.
      q_b3_low = map_dbl(Est_HMM_boot, function(x) {
        base <- mutate(x,
                       eps =   (Real_boot - V_HMM_small_est_boot) / Size_boot)
       quant_max( alpha, base)
      }),
  
      q_aldemi_high = map_dbl(Est_HMM_boot, function(x) {
        base <- mutate(x,
                       eps =   (V_HMM_oracle_boot_aldemi - V_HMM_est_boot_aldemi) / Size_boot)
        quant_max(1 - alpha / 2, base)
      }),
      q_aldemi_low = map_dbl(Est_HMM_boot, function(x) {
        mutate(
          x,
          eps =   (
            V_HMM_small_oracle_boot_aldemi - V_HMM_small_est_boot_aldemi
          ) / Size_boot
        ) %>%
          quant_max(alpha / 2, .)
      }),
      q_aldelta_high = map_dbl(Est_HMM_boot, function(x) {
        base <- mutate(x,
                       eps =   (V_HMM_oracle_boot_aldelta - V_HMM_est_boot_aldelta) / Size_boot)
        quant_max(1 - alpha *(1-delt), base)
      }),
      q_aldelta_low = map_dbl(Est_HMM_boot, function(x) {
        mutate(
          x,
          eps =   (
            V_HMM_small_oracle_boot_aldelta - V_HMM_small_est_boot_aldelta
          ) / Size_boot
        ) %>%
          quant_max(alpha * (1-delt), .)
      }),
      q_al1_moins_delta_high = map_dbl(Est_HMM_boot, function(x) {
        base <- mutate(x,
                       eps =   (V_HMM_oracle_boot_al1_moins_delta - V_HMM_est_boot_al1_moins_delta) / Size_boot)
        quant_max(1 - alpha *( delt), base)
      }),
      q_al1_moins_delta_low = map_dbl(Est_HMM_boot, function(x) {
        mutate(
          x,
          eps =   (
            V_HMM_small_oracle_boot_al1_moins_delta - V_HMM_small_est_boot_al1_moins_delta
          ) / Size_boot
        ) %>%
          quant_max(alpha * (delt), .)
      }),
      q_samesel = map2_dbl(Est_HMM_boot, V_HMM_est_aldemi, function(x, y) {
        mutate(x, eps =  (y -  V_HMM_est_boot_aldemi_samesel)/Size) %>%
          quant_max(1 - alpha / 2, .)
          # mutate(eps =  map_dbl(eps, ~ min(max(.x, 0), Size))) %>%
          # pull(eps) %>%
          # quantile(probs = 1 - alpha / 2,
          #          na.rm = TRUE,
          #          type = 3)
        
      }),
      q_samesel_small = map2_dbl(Est_HMM_boot, V_HMM_small_est_aldemi, function(x, y) {
        mutate(x, eps =  (y -  V_HMM_small_est_boot_aldemi_samesel) / Size) %>%
          quant_max(alpha / 2, .)
        
      }),
      q_samesel_delta = map2_dbl(Est_HMM_boot, V_HMM_est_aldelta, function(x, y) {
        mutate(x, eps =  (y -  V_HMM_est_boot_aldelta_samesel)/Size) %>%
          quant_max(1 - alpha *(1- delt), .)
        # mutate(eps =  map_dbl(eps, ~ min(max(.x, 0), Size))) %>%
        # pull(eps) %>%
        # quantile(probs = 1 - alpha / 2,
        #          na.rm = TRUE,
        #          type = 3)
        
      }),
      q_samesel_small_delta = map2_dbl(Est_HMM_boot, V_HMM_small_est_aldelta, function(x, y) {
        mutate(x, eps =  (y -  V_HMM_small_est_boot_aldelta_samesel) / Size) %>%
          quant_max(alpha *(1- delt), .)
        
      }),
      q_samesel1_moins_delta = map2_dbl(Est_HMM_boot, V_HMM_est_al1_moins_delta, function(x, y) {
        mutate(x, eps =  (y -  V_HMM_est_boot_al1_moins_delta_samesel)) %>%
          quant_max(1 - alpha *( delt), .)
      }),
     
       q_samesel_small1_moins_delta = map2_dbl(Est_HMM_boot, V_HMM_small_est_al1_moins_delta, function(x, y) {
        mutate(x, eps =  (y -  V_HMM_small_est_boot_al1_moins_delta_samesel) / Size) %>%
          quant_max(alpha * delt, .)
        
      }),
      V_HMM_boot_samesel = V_HMM_est_aldemi + q_samesel * Size,
      V_HMM_small_boot_samesel = V_HMM_small_est_aldemi + q_samesel_small *
        Size,
      V_HMM_boot_samesel_delta = V_HMM_est_aldelta + q_samesel_delta * Size ,
      V_HMM_small_boot_samesel_delta = V_HMM_small_est_aldelta + q_samesel_small_delta *
        Size,
      V_HMM_boot3 = V_HMM_est + q_b3_high * Size,
      V_HMM_small_boot3 = V_HMM_small_est + q_b3_low * Size,
      V_HMM_boot_q1_moins_delta = V_HMM_est_al1_moins_delta + q_al1_moins_delta_high * Size,
      V_HMM_small_boot_q1_moins_delta = V_HMM_small_est_al1_moins_delta + q_al1_moins_delta_low * Size,
      
      V_HMM_boot_qdelta = V_HMM_est_aldelta + q_aldelta_high * Size,
      V_HMM_small_boot_qdelta = V_HMM_small_est_aldelta + q_aldelta_low * Size,
      V_HMM_boot_qdemi = V_HMM_est_aldemi + q_aldemi_high * Size,
      V_HMM_small_boot_qdemi = V_HMM_small_est_aldemi + q_aldemi_low * Size,
      V_HMM_boot_naif = map_dbl(
        Est_HMM_boot,
        ~ mutate(., eps = Real_boot / Size_boot) %>%
          quant_max((1-alpha), .)
      ),
      V_HMM_boot_naif = V_HMM_boot_naif * Size,
      V_HMM_small_boot_naif = map_dbl(
        Est_HMM_boot,
        ~  mutate(., eps = Real_boot / Size_boot) %>%
          quant_max((alpha), .)
      ),
      V_HMM_small_boot_naif = V_HMM_small_boot_naif * Size,
      FDR_boot_naif = map_dbl(Est_HMM_boot, ~ pull(. , Real_boot) %>%
                                mean()),
      Nom = case_when(
        Nom == "pval_tresh" ~  "bgroup('{',p < 0.05,'}')",
        Nom == "pvalm_sel" ~   "bgroup('{',p < 0.05,'}')~intersect()~bgroup('{','1, ..., 200' ,'}')",
        Nom == "sel_viter_est" ~ "Viterbi",
        Nom == "sel_viter_min_size" ~ "Viterbi + 4",
        Nom == "pvalm_tresh_H0" ~  "H0~intersect()~ bgroup('{',p < 0.05 ,'}')",
        Nom == "sel_viter_est_H0" ~ "H0~intersect()~Viterbi",
        Nom == "H1" ~  "H1",
        Nom == "pval_lfdr" ~  "bgroup('{',p < th,'}')",
        Nom == "lfdr_tresh" ~  "SC(0.05)",
        Nom == "sel_idpval" ~  "SC(FDR[p])",
        Nom == "lfdr_tresh_viter" ~  "SC(FDR[v])",
        Nom == "block" ~ "bgroup('{','1, ..., 200' ,'}')"
      ),
      Nom = factor(
        Nom,
        levels = c(
          "bgroup('{',p < 0.05,'}')",
          "bgroup('{',p < 0.05,'}')~intersect()~bgroup('{','1, ..., 200' ,'}')",
          "SC(FDR[p])",
          "SC(0.05)",
          "bgroup('{',p < th,'}')",
          "Viterbi",
          "SC(FDR[v])",
          "Viterbi + 4",
          "bgroup('{','1, ..., 200' ,'}')",
          "H1",
          "H0~intersect()~ bgroup('{',p < 0.05 ,'}')",
          "H0~intersect()~Viterbi"
        )
      )
    ) %>%
    mutate_at(.vars = vars(starts_with("V_")), ~ case_when(. > Size ~  Size,
                                                           TRUE ~ .))
  return(Res_all)
}


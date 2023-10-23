#' Title
#'
#' @param Final output of sanssouci.hmm
#' @param sim logical whether it come from a simulation (if TRUE de TRUE FDP is represented by a cross)
#'
#' @return
#' @export
#'
#' @examples
plot_IC  <- function(Final, col1 = c("orange","lightblue","grey30"), sim = FALSE) {
  if(!sim){
  IC <- Final %>%
    filter(!is.na(Nom) & as.numeric(Nom)!=8) %>%
    mutate_at(vars(starts_with("V")),~./Size) %>% 
    mutate(
      FDR_est = FDR_est/Size,
      V_HMM_boot_qdemi = case_when(
        is.na(V_HMM_boot_qdemi)~ V_HMM_boot_samesel,
        TRUE ~ V_HMM_boot_qdemi
      ),
      V_HMM_small_boot_qdemi = case_when(
        is.na(V_HMM_small_boot_qdemi)~ V_HMM_small_boot_samesel,
        TRUE ~ V_HMM_small_boot_qdemi
      )
    ) %>%
    select(Nom,
           # starts_with("V_simes"), V_small_simes,
           
           V_HMM_small_boot3, V_HMM_boot3,
           V_HMM_small_boot_samesel, V_HMM_boot_samesel,
           # V_HMM_small_boot_qdemi, V_HMM_boot_qdemi,
           V_HMM_small_est, V_HMM_est,
           FDR_est) %>%
    gather(-Nom, -FDR_est,
           # -V_HMM_boot_qdemi,
           -V_HMM_est, -V_HMM_boot3,-V_HMM_boot_samesel,
           # -V_simes,
           value ="high", key ="Type1") %>%
    gather(-Nom, -FDR_est, -high, -Type1,
           value ="small", key ="Type2") %>%
    mutate(Type1 = gsub("small_","", Type1),
           Type1 = factor( Type1,
                           levels =
                             c("V_HMM_boot3","V_HMM_boot_samesel",
                               "V_HMM_est"))
           ) %>%
    # mutate(Type2 = gsub("V_", "", Type2)) %>%
    filter(Type1 == Type2) %>%
    ggplot(aes(x = Nom, ymin = small,
               ymax = high, color = Type1, size = Type1)) +
    geom_errorbar(alpha = 0.6) + coord_flip()+
    geom_point(aes(y = FDR_est))+
    scale_x_discrete(labels = function(l)parse(text=l)) +
    labs(x ="", y = "FDP") +
    scale_color_manual(values =  col1
                       ,labels =c("Boot3","Boot2","PI")) +
    scale_size_manual(values = c(3,2.5,1.2),
                      labels =c("Boot3","Boot2","PI"))+
    theme_bw() +
    theme(legend.position = "bottom") +
    guides(color=guide_legend(nrow=1,byrow=TRUE)) 
} else{
  IC <- Final %>%
    filter(!is.na(Nom) & as.numeric(Nom)!=8) %>%
    mutate_at(vars(starts_with("V")),~./Size) %>% 
    mutate(
      FDR_est = FDR_est/Size,
      V_HMM_boot_qdemi = case_when(
        is.na(V_HMM_boot_qdemi)~ V_HMM_boot_samesel,
        TRUE ~ V_HMM_boot_qdemi
      ),
      V_HMM_small_boot_qdemi = case_when(
        is.na(V_HMM_small_boot_qdemi)~ V_HMM_small_boot_samesel,
        TRUE ~ V_HMM_small_boot_qdemi
      )
    ) %>%
    select(Nom,
           # starts_with("V_simes"), V_small_simes,
           
           V_HMM_small_boot3, V_HMM_boot3,
           V_HMM_small_boot_samesel, V_HMM_boot_samesel,
           # V_HMM_small_boot_qdemi, V_HMM_boot_qdemi,
           V_HMM_small_est, V_HMM_est,
           FDR_est, V_real) %>%
    gather(-Nom, -FDR_est, -V_real,
           # -V_HMM_boot_qdemi,
           -V_HMM_est, -V_HMM_boot3,-V_HMM_boot_samesel,
           # -V_simes,
           value ="high", key ="Type1") %>%
    gather(-Nom, -FDR_est,-V_real, -high, -Type1,
           value ="small", key ="Type2") %>%
    mutate(Type1 = gsub("small_","", Type1),
           Type1 = factor( Type1,
                           levels =
                             c("V_HMM_boot3","V_HMM_boot_samesel",
                               "V_HMM_est"))
    ) %>%
    # mutate(Type2 = gsub("V_", "", Type2)) %>%
    filter(Type1 == Type2) %>%
    ggplot(aes(x = Nom, ymin = small,
               ymax = high, color = Type1, size = Type1)) +
    geom_errorbar(alpha = 0.6) + coord_flip()+
    geom_point(aes(y = FDR_est))+geom_point(aes(y = V_real), shape = 4, size = 3)+
    scale_x_discrete(labels = function(l)parse(text=l)) +
    labs(x ="", y = "FDP") +
    scale_color_manual(values =  col1
                       ,labels =c("Boot3","Boot2","PI")) +
    scale_size_manual(values = c(3,2.5,1.2),
                      labels =c("Boot3","Boot2","PI"))+
    theme_bw() +
    theme(legend.position = "bottom") +
    guides(color=guide_legend(nrow=1,byrow=TRUE)) 
  
}
IC
}






#' Title
#'
#' @param Final 
#'
#' @return
#' @export
#'
#' @examples
plot_IC  <- function(Final, col1 = c("orange","lightblue","grey30")) {
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
  IC
}



plot_IC_simu  <- function(Final) {
  browser()
  Tb <- Final %>%
    mutate(FDR_qdemi = FDR_est) %>%
    nest(
      high = c(V_HMM_or, V_HMM_boot_samesel, V_HMM_boot3),
      law = c(V_HMM_small_or, V_HMM_small_boot_samesel,
              V_HMM_small_boot3),
      FDR = c(FDR_or, FDR_est, FDR_qdemi)
    ) %>%
    mutate(
      high = map(high, ~ gather(., key = "high_type", value = "high")),
      law = map(law, ~ gather(., key = "low_type", value = "low")),
      FDR = map(FDR, ~ gather(., key = "FDR_type", value = "FDR"))
    ) %>%
    unnest(cols = c(high, law, FDR)) %>%
    mutate(
      FDR_type  = factor(FDR_type, levels = c("FDR_or", "FDR_est", "FDR_qdemi")),
      high = high / Size,
      low = low / Size,
      FDR = FDR / Size,
      Real = V_real / Size,
      Runs = paste0("runs", simu)
    )
  table(Tb$low_type, Tb$FDR_type)
  Real_FDR <- Tb %>% select(FDR, FDR_type, Nom,  Real,Runs) %>%
    gather(key = "type", value = "value",-FDR_type,-Nom,-Runs) %>%
    mutate(Col = case_when(type == "FDR" ~ "blue",
                           TRUE ~ "red"))
  # %>%
  #   gather(key = "type", value ="value", - FDR_type, -Nom, -simu) %>%
  #   mutate(Col = case_when( type == "FDR" ~ "blue",
  #                           TRUE ~ "red"))
  
  med_FDR <- Tb %>% select(FDR, FDR_type, Nom, Runs)
  real <- Tb %>% select(FDR_type, Nom, Runs, Real)
  # %>%
  #   gather(key = "type", value ="value", - FDR_type, -Nom, -simu)
  Tb  %>% 
    ggplot(aes(x = FDR_type,
               group = Runs))  +
    scale_x_discrete(
      labels = function(l)
        parse(text = l)
    ) +
    geom_errorbar(aes(
      ymin = low,
      ymax = high,
      color = FDR_type
    ),
    position = position_dodge(width = 1))  +
    scale_color_manual(values = col1, labels =
                         c("Oracle", "boot2", "boot3")) +
    # new_scale_color() +
    
    geom_point(
      data = med_FDR,
      aes(y = FDR, color = FDR_type),
      position = position_dodge(width = 1),
      shape = 2,
      size = 2
    ) +
    geom_point(
      data = real,
      aes(y = Real, color = FDR_type),
      position = position_dodge(width = 1),
      shape = 4,
      size = 2
    ) +
    geom_point(data = Real_FDR, y = NA , aes(shape = type)) +
    scale_shape_manual(values = c(2, 4), labels = c("FDR estimate", "Actual FDP")) +
    facet_grid( ~ Nom, labeller = label_parsed) + coord_flip() + 
    theme(axis.text.y = element_blank())+labs(x="")
  
}


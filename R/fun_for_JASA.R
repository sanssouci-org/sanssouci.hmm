## This is the function to get the graphe of the paper.


## First function : get the bootstrap result
# and the differents bounds


graph_diff_high_FDR <-
  function(Res,
           order,
           type = "high",
           num_1 = 0.5,
           num_2 = 0.5,
           row_wrap = 1) {
    f_type <- function(x) {
      sum(x < 0, na.rm = TRUE) / length(x)
    }
    if (type == "high") {

      Diff <-
        Res %>%  select(all_of(c(order, "Nom", "V_real", "Size"))) %>%
        mutate(real = V_real) %>%
        mutate_at(.vars = vars(starts_with("V_")), ~ as.numeric(. - real)) %>%
        select(-V_real, -real)
      
      
    }
    if (type == "low") {
        Diff <-
        Res %>%  select(all_of(c(order, "Nom", "V_real", "Size"))) %>%
        mutate(real = V_real) %>%
        mutate_at(.vars = vars(starts_with("V_")), ~ as.numeric( real - .)) %>%
        select(-V_real, -real)
      
      
    }
    Data <- Diff %>%
      mutate_at(all_of(order),  ~ . / Size) %>%
      select(all_of(c(order, "Nom"))) %>%
      gather(-Nom, key = "key",
             value = "value") 
    N1 <- Data %>%
      group_by(Nom) %>%
      summarise(y=max(value, na.rm=TRUE)) 
    Prop <-  Diff %>%
      group_by(Nom) %>%
      summarise_at(vars(starts_with("V_")),  f_type) %>%
      gather(-Nom, key = "key",
             value = "labs") %>%
      mutate(
        labs = round(labs * 100),
        labs = paste(labs, "%")
        # y = case_when(as.numeric(Nom) %in% c(1:3, 5, 7, 9) ~ num_1,
        #               TRUE ~ num_1)
      ) %>% 
      left_join(N1, by = "Nom") %>% filter(labs != "NA %")
  
    p <- Data %>%
      mutate(key = factor(key, levels = order)) %>%
      ggplot(aes(x = key, y = value, color = key)) +
      
      geom_jitter(alpha = 0.5) +
      # expression(bgroup("|",H[0],"|") -  V^{HMM})
      coord_flip() +
      labs(y = "",
           x = "") +
      theme(
        legend.title = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "right"
      ) +
      geom_hline(yintercept = 0, linetype ="dotted") +
      facet_wrap( ~ Nom,
                  scale = "free",
                  labeller = label_parsed,
                  nrow = row_wrap) +
      geom_label(data = Prop,
                 aes(label = labs, y = y),
                 show.legend = F)  +
      scale_x_discrete(limits = rev(order))
    p
  }

## Graphe pour la puissance 
graphe_puiss_mean <- function(Res, order){
  f_type <- function(x) {
    sum(x < 0) / length(x)
  }
  Diff <-
    Res %>%  select(all_of(c(order, "Nom", "V_real","Size"))) %>%
    mutate(real = V_real) %>%
    mutate_at(.vars = vars(starts_with("V_")), ~ (Size-.) / (Size -real)) %>%
    select( -real, -Size)
  
  p <- Diff %>%
    select(all_of(c(order, "Nom"))) %>%
    gather(-Nom, key = "key",
           value = "value") %>%
    filter(value <=1) %>% 
    mutate(key = factor(key, levels = order)) %>%
    ggplot(aes(x = key, y = value, color = key)) +
    geom_violin() +
    stat_summary() +
    # expression(bgroup("|",H[0],"|") -  V^{HMM})
    coord_flip() +
    labs(y = "",
         x = "") +
    theme(
      legend.title = element_blank(),
      axis.text.y = element_blank(),
      legend.position = "right"
    ) +
    facet_wrap( ~ Nom,
                scale = "free",
                labeller = label_parsed,
                nrow = 1) +
    scale_x_discrete(limits = rev(order))
  p
}

graph_diff_high_FDR_facet <-
  function(Res,
           order,
           type = "high",
           num_1 = 0.5,
           num_2 = 0.5,
           row_wrap = 1,
           var_comp ) {
    if (type == "high") {
      f_type <- function(x) {
        sum(x < 0) / length(x)
      }
    }
    if (type == "low") {
      f_type <- function(x) {
        sum(x > 0) / length(x)
      }
    }
    Diff <-
      Res %>%  select(all_of(c(order, "Nom", "V_real", var_comp, "Size"))) %>%
      mutate(real = V_real, Var_comp = !!sym(var_comp)) %>%
      mutate_at(.vars = vars(starts_with("V_")), ~ as.numeric(. - real)) %>%
      select(-V_real, -real)
    
    
    Prop <-  Diff %>%
      group_by(Nom, Var_comp) %>%
      summarise_at(vars(starts_with("V_")),  f_type) %>%
      gather(-Nom,-Var_comp, key = "key",
             value = "labs") %>%
      mutate(
        labs = round(labs * 100),
        labs = paste(labs, "%"),
        y = case_when(as.numeric(Nom) %in% c(1:3, 5, 7, 9) ~ num_1,
                      TRUE ~ num_2)
      )
    p <- Diff %>%
      mutate_at(all_of(order),  ~ . / Size) %>%
      select(all_of(c(order, "Nom","Var_comp"))) %>%
      gather(-Nom,-Var_comp, key = "key",
             value = "value") %>%
      mutate(key = factor(key, levels = order)) %>%
      ggplot(aes(x = key, y = value, color = key)) +
      
      geom_jitter(alpha = 0.5) +
      # expression(bgroup("|",H[0],"|") -  V^{HMM})
      coord_flip() +
      labs(y = "",
           x = "") +
      theme(
        legend.title = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "right"
      ) +
      geom_hline(yintercept = 0, linetype ="dotted") +
      facet_grid(Var_comp ~ Nom,
                 scale = "free",
                 labeller = label_parsed) +
      geom_label(data = Prop,
                 aes(label = labs, y = y),
                 show.legend = F)  +
      scale_x_discrete(limits = rev(order)) 
      p
  }
graphe_puiss_mean_facet <- function(Res, order, var_comp){
  f_type <- function(x) {
    sum(x < 0) / length(x)
  }
  Diff <-
    Res %>%  select(all_of(c(order, "Nom", "V_real","Size",var_comp))) %>%
    mutate(real = V_real,  Var_comp = !!sym(var_comp)) %>%
    mutate_at(.vars = vars(starts_with("V_")), ~ (Size-.) / (Size -real)) %>%
    select( -real, -Size)
  
  p <- Diff %>%
    select(all_of(c(order, "Nom","Var_comp"))) %>%
    gather(-Nom,-Var_comp, key = "key",
           value = "value") %>%
    filter(value <=1) %>% 
    mutate(key = factor(key, levels = order)) %>%
    ggplot(aes(x = key, y = value, color = key)) +
    geom_violin() +
    stat_summary() +
    # expression(bgroup("|",H[0],"|") -  V^{HMM})
    coord_flip() +
    labs(y = "",
         x = "") +
    theme(
      legend.title = element_blank(),
      axis.text.y = element_blank(),
      legend.position = "right"
    ) +
    facet_grid(Var_comp ~ Nom,
                scale = "free",
                labeller = label_parsed) +
    scale_x_discrete(limits = rev(order))
  p
}




bx_plot_prop <- function(Res_all, num_sim = 5, col1) {
  Tb <- Res_all %>%
    filter(simu %in% 1:num_sim) %>%
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


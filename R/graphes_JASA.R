## graphe for paper 

rm(list=ls())
library(tidyverse)
library(sansSouci)
library(ggsci)
source('fun_for_JASA.R')
source("graphes.R")
library(extrafont) 
library(viridis)
theme_set(theme_bw() +
            theme(strip.background = element_rect(fill = "white"),
                  text = element_text(face="bold", family="LM Roman 10"),
            ))

require(ggsci)
get_res <- function(name){
  # l <- load(paste0(name, ".RData"))
  # Result <- get(l)
  # res<- Result %>% 
  #   # get_results()  
  #   mise_en_forme()
  # save(res, file = paste0("res_",name, ".Rdata"))
  # res
}





plot_sim <- function(name, type="sel"){
  l <- load(paste0("res_",name, ".Rdata"))
  res <- get(l)
  res <- res %>% 
    filter(as.numeric(Nom)!=8)
  # col1 <- viridis_pal()(9)
  cols <- pal_nejm()(8)
  col1 <- c(cols[1:2],c("chocolate1","chocolate2","chocolate3"),
            cols[c(4,7,5)])
  sel <- c(1,2,7,8)
  lab_high <-  c(  
                   "Naive",
                   "Oracle",
                   "Plug-in" , 
                   "Boot1 0.9",
                   "Boot1 0.5",
                   "Boot1 0.1",
                   "Boot2",
                   "Boot3")
  ord_high <- c( 
                 "V_HMM_boot_naif",
                 "V_HMM_or",
                 "V_HMM_est",
                 "V_HMM_boot_qdelta",
                 "V_HMM_boot_qdemi",
                 "V_HMM_boot_q1_moins_delta",
                 "V_HMM_boot_samesel", 
                 "V_HMM_boot3")
  if(type!="sel"){
    p1 <- res  %>% 
      select( - V_HMM_est_aldemi ) %>% 
      graph_diff_high_FDR(order = rev(ord_high), row_wrap =1, 
                          num_1 = 0.1, num_2 =0.1) + 
      scale_color_manual(values = col1, labels = rev(lab_high)) +
      theme(legend.text = element_text(size=14))+
      guides(color=guide_legend(override.aes = list(size=3, alpha =1)))
  } else{
    p1 <- res  %>% 
      select( - V_HMM_est_aldemi ) %>% 
      graph_diff_high_FDR(order = rev(ord_high[sel]), row_wrap =1,
                          num_1 = 0.1, num_2 =0.1) + 
      scale_color_manual(values = col1[sel], labels = rev(lab_high[sel])) +
      theme(legend.text = element_text(size=14))+
      guides(color=guide_legend(override.aes = list(size=3, alpha =1)))
    
  }
  ggsave(p1, file =paste0(name,"_diff.pdf"), device =cairo_pdf, width= 16, height =4)
  
  
  ## graphe de puissance : 
  
  p2 <- res  %>% 
    select( - V_HMM_est_aldemi ) %>% 
    graphe_puiss_mean(order = rev(ord_high[sel])) + 
    scale_color_manual(values = col1[sel], labels = rev(lab_high[sel])) +
    theme(legend.text = element_text(size=14))
  ggsave(p2, file =paste0(name,"puiss_mean.pdf"), device =cairo_pdf, width= 16, height =4)
  
  lab_low <-  c(  "Naive",
                  "Oracle",
                  "Boot2",
                  "Boot3")
  ord_low <- c(
    "V_HMM_small_boot_naif",
    "V_HMM_small_or",
    "V_HMM_small_boot_samesel", 
    "V_HMM_small_boot3")
  
  p1 <- res  %>% 
    select( - V_HMM_est_aldemi ) %>% 
    # mutate(Size = as.character(Size)) %>% 
    # mutate_if(is.numeric, ~./as.numeric(Size)) %>% 
    graph_diff_high_FDR(order = rev(ord_low),
                        row_wrap =1, num_1 = 0.07, 
                        num_2 = 0.07, type ="low") + 
    # scale_color_manual(labels = rev(lab_low[-1]), values = col1) +
    scale_color_manual(values = col1[sel],labels = rev(lab_low)) +
    theme(legend.text = element_text(size=14))+
    guides(color=guide_legend(override.aes = list(size=3, alpha =1)))
  ggsave(p1, file =paste0(name, "_lower_diff",".pdf"), device =cairo_pdf, width= 16, height =4)
  
  IC <-res  %>% 
    select( - V_HMM_est_aldemi ) %>% 
    bx_plot_prop(num_sim =3, col1 = rev(col1[sel[-4]]))+
    scale_color_manual(values = rev(col1[sel[-4]]),labels = lab_low[-1]) +
    theme(legend.text = element_text(size=14),legend.title = element_blank())+
    guides(color=guide_legend(override.aes = list(size=3, alpha =1)))+
    scale_shape_manual(values = c(2, 4), labels = c("FDR estimate", "Actual FDP")) 
  ggsave(IC, file =paste0(name, "_IC",".pdf"), device =cairo_pdf, width= 16, height =4)
  
  # return(res)
}



 Res <- get_res("res_simu_delta")

plot_sim("res_simu_delta", type ="all")

## avec H0 knowledge 
 Res <- get_res("res_simu_delta_h0")

plot_sim("res_simu_delta_h0")

## unknown f0 : 


 res_uf0 <- get_res("delta_highdet_m3200_uf0")

plot_sim("delta_highdet_m3200_uf0")
# 
  # res_ff0 <- get_res("delta_highdet_uf0_start_falsef0_m3200")
# 
# plot_sim("delta_highdet_uf0_start_falsef0_m3200")

  Res <- get_res("delta_highdet_uf0_start_locfdr_m3200")

plot_sim("delta_highdet_uf0_start_locfdr_m3200")
##
## graphe realiste


## graphe avec les determinants : 

  # res_realiste<- get_res("Realist_frac1")
# 
# plot_sim("Realist_frac1")
# 
  Res <- get_res("delta_diffdet_m3200")
cols <- pal_nejm()(8)
col1 <- c(cols[1:2],c("chocolate1","chocolate2","chocolate3"),
          cols[c(4,7,5)])

lab_high_tot <-  c(  
  "Naive",
  "Oracle",
  "Plug-in" , 
  "Boot1 0.9",
  "Boot1 0.5",
  "Boot1 0.1",
  "Boot2",
  "Boot3")
ord_high_tot <- c( 
  "V_HMM_boot_naif",
  "V_HMM_or",
  "V_HMM_est",
  "V_HMM_boot_qdelta",
  "V_HMM_boot_qdemi",
  "V_HMM_boot_q1_moins_delta",
  "V_HMM_boot_samesel", 
  "V_HMM_boot3")
l <- load(paste0("res_delta_diffdet_m3200.Rdata"))
res <- get(l)
sel <- c(1,2,7,8)
res  %>% 
  mutate(Det_A = map_dbl(A, ~round(det(.),2))) %>% 
  filter(Det_A==0) %>% 
  mutate(Diff= V_real -V_HMM_boot_qdemi) %>% 
  pull(Diff) %>% quantile(probs = seq(0,1,0.1),na.rm=TRUE)

lab_high <-  c("Simes",lab_high_tot[sel])
ord_high <- c("V_simes",
              ord_high_tot[sel])
p1 <- res  %>% 
  mutate(Det_A = map_chr(A, ~paste("Det(A) ==",round(det(.),2)))) %>% 
  select( - V_HMM_est_aldemi ) %>% 
  graph_diff_high_FDR_facet(order = rev(ord_high), 
                            row_wrap =1, var_comp ="Det_A") + 
  scale_color_manual(labels = rev(lab_high), values = c(col1[sel],"grey50")) +
  theme(legend.text = element_text(size=14))+
  guides(color=guide_legend(override.aes = list(size=3, alpha =1)))
ggsave(p1, file =paste0("delta_diffdet_m3200_diff.pdf"),
       device =cairo_pdf, width= 16, height =12)
# return(res)


  res_realiste<- get_res("Realist_frac1")

  res_realiste_otherfrac<- get_res("Realist_other_frac_delta")

load("res_Realist_other_frac_delta.Rdata")
res_realiste_otherfrac<-res
load("res_Realist_frac1.Rdata")
sel <- c(1,7,8)
sel_col <-c(1,2,8)
res_reali <- rbind(res, res_realiste_otherfrac)
lab_high <-  c(  "Simes",lab_high_tot[sel])
ord_high <- c("V_simes",ord_high_tot[sel])
p1 <- res_reali  %>% 
  filter(!(as.numeric(Nom)%in% c(2,8))) %>% 
  mutate(Det_A = map_chr(tumorFraction, ~paste("Tf ==",.)), 
         Det_A = factor(Det_A, 
                        levels = c("Tf == 1", "Tf == 0.7","Tf == 0.5","Tf == 0.3") )) %>% 
  select( - V_HMM_est_aldemi ) %>% 
  graph_diff_high_FDR_facet(order = rev(ord_high),
                            row_wrap =1, var_comp ="Det_A") + 
  scale_color_manual(values = c(col1[sel_col], "grey50"), labels = rev(lab_high)) +
  theme(legend.text = element_text(size=14))+
  guides(color=guide_legend(override.aes = list(size=3, alpha =1))) +
  scale_y_continuous(breaks = c(0,0.2,0.5,0.7,1))
ggsave(p1, file =paste0("reali_m10000_diff.pdf"),
       device =cairo_pdf, width= 10, height =6)

## puissance 
p1 <- res_reali  %>% 
  filter(!(as.numeric(Nom)%in% c(2,8))) %>% 
  mutate(Det_A = map_chr(tumorFraction, ~paste("Tf ==",.)), 
         Det_A = factor(Det_A, 
                        levels = c("Tf == 1", "Tf == 0.7","Tf == 0.5","Tf == 0.3") )) %>% 
  select( - V_HMM_est_aldemi ) %>% 
  graphe_puiss_mean_facet(order = rev(ord_high), var_comp ="Det_A") + 
  scale_color_manual(values = c(col1[sel_col], "grey50"), labels = rev(lab_high)) +
  theme(legend.text = element_text(size=14))+
  guides(color=guide_legend(override.aes = list(size=1, alpha =1))) +
  scale_y_continuous(breaks = c(0,0.2,0.5,0.7,1))
ggsave(p1, file =paste0("reali_m10000_puiss.pdf"),
       device =cairo_pdf, width= 10, height =6)


## Ajout BNR et Bram : 
load("~/posthoc/HMM/articleMarie/JASA/Simu_papier_marie/Res_simu_other.RData")

name <- "res_simu_delta"

l <- load(paste0("res_",name, ".Rdata"))
res <- get(l) %>% 
  filter(as.numeric(Nom)!=8)
res2 <- Res_simu_other %>% mutate( Nom = case_when(
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
)) %>% 
  filter(as.numeric(Nom)!=8)
res <- res %>% mutate(Nom = as.character(Nom)) %>%  left_join(res2) 
# col1 <- viridis_pal()(9)
cols <- pal_nejm()(8)
col1 <- c(rep("grey",3),"grey20", cols[1:2],c("chocolate1","chocolate2","chocolate3"),
          cols[c(4,7,5)])
sel <- c(1,2,7,8)
lab_high <-  c(  
  "Naive",
  "Oracle",
  "Plug-in" , 
  "Boot1 0.9",
  "Boot1 0.5",
  "Boot1 0.1",
  "Boot2",
  "Boot3",
  "Katsevitch, Ramdas ",
  "BNR m",
  "BNR 2m1","Simes")
ord_high <- c( 
  "V_HMM_boot_naif",
  "V_HMM_or",
  "V_HMM_est",
  "V_HMM_boot_qdelta",
  "V_HMM_boot_qdemi",
  "V_HMM_boot_q1_moins_delta",
  "V_HMM_boot_samesel", 
  "V_HMM_boot3",
  "V_bram",
  "V_bnr_m",
  "V_bnr_2m1", "V_simes")
  p1 <- res  %>% 
    select( - V_HMM_est_aldemi ) %>% 
    graph_diff_high_FDR(order = rev(ord_high), row_wrap =1, 
                        num_1 = 0.1, num_2 =0.1) + 
    scale_color_manual(values = col1, labels = rev(lab_high)) +
    theme(legend.text = element_text(size=14))+
    guides(color=guide_legend(override.aes = list(size=3, alpha =1)))

  p1
ggsave(p1, file =paste0("compare_BNR_ramdas_diff_all.pdf"), device =cairo_pdf, width= 16, height =4)

## graphe pour introduction :
library(ggthemes)
col_intro <- c(  few_pal()(4),"grey20",col1[c(5,10,11)])
lab_high <-  c(  
  "Oracle",
  "Plug-in" , 
  "Boot3",
  "KR",
  "BNR m",
  "DBNR",
  "DBNR_tree",
  "Simes")
ord_high <- c( 
  "V_HMM_or",
  "V_HMM_est",
  "V_HMM_boot3",
  "V_bram",
  "V_bnr_m",
  "V_DBNR_part",
  "V_DBNR_tree",
  "V_simes")
p1 <- res  %>% 
  select( - V_HMM_est_aldemi ) %>% 
  graph_diff_high_FDR(order = rev(ord_high), row_wrap =1, 
                      num_1 = 0.1, num_2 =0.1) + 
  scale_color_manual(values = col_intro, labels = rev(lab_high)) +
  theme(legend.text = element_text(size=14))+
  guides(color=guide_legend(override.aes = list(size=3, alpha =1)))

p1
ggsave(p1, file =paste0("compare_intro.jpeg"), 
        width= 16, height =4)


col_intro <- c(  few_pal()(4)[3:4],"grey20",col1[c(5,11)])
col_intro <- c("mediumblue", "slateblue2", "orchid3",col1[c(5,11)])
col_intro <- c("paleturquoise3", "palevioletred", "grey30",col1[c(5,10,11)])

lab_high <-  c(  
  "Oracle",
  "Plug-in",
  "Boot3",
  "KR",
  "BNR","Simes")
ord_high <- c( 
  "V_HMM_or",
  "V_HMM_est",
  "V_HMM_boot3",
  "V_bram",
  "V_bnr_m",
  "V_simes")
p1 <- res  %>% 
  select( - V_HMM_est_aldemi ) %>% 
  graph_diff_high_FDR(order = rev(ord_high), row_wrap =1, 
                      num_1 = 0.1, num_2 =0.1) + 
  scale_color_manual(values = col_intro, labels = rev(lab_high)) +
  theme(text = element_text(size=14))+
  guides(color=guide_legend(override.aes = list(size=3, alpha =1)))

p1
ggsave(p1, file =paste0("compare_intro_ssdbnr.pdf"),device =cairo_pdf, 
       width= 13, height =3)

## en vrai je pense pas mettre simes et BNR 

cols <- pal_nejm()(8)
col1 <- c("grey20", cols[1:2],c("chocolate1","chocolate2","chocolate3"),
          cols[c(4,7,5)])
sel <- c(1,2,7,8)
lab_high <-  c(  
  "Naive",
  "Oracle",
  "Plug-in" , 
  "Boot1 0.9",
  "Boot1 0.5",
  "Boot1 0.1",
  "Boot2",
  "Boot3",
  "KR")
ord_high <- c( 
  "V_HMM_boot_naif",
  "V_HMM_or",
  "V_HMM_est",
  "V_HMM_boot_qdelta",
  "V_HMM_boot_qdemi",
  "V_HMM_boot_q1_moins_delta",
  "V_HMM_boot_samesel", 
  "V_HMM_boot3",
  "V_bram")
p1 <- res  %>% 
  select( - V_HMM_est_aldemi ) %>% 
  graph_diff_high_FDR(order = rev(ord_high), row_wrap =1, 
                      num_1 = 0.1, num_2 =0.1) + 
  scale_color_manual(values = col1, labels = rev(lab_high)) +
  theme(legend.text = element_text(size=14))+
  guides(color=guide_legend(override.aes = list(size=3, alpha =1)))
p1
ggsave(p1, file =paste0("compare_BNR_ramdas_diff.pdf"), device =cairo_pdf, width= 16, height =4)
sel <- c(1,2,7,8,9)
p2 <- res  %>% 
  select( - V_HMM_est_aldemi ) %>% 
  graphe_puiss_mean(order = rev(ord_high[sel])) + 
  scale_color_manual(values = col1[sel], labels = rev(lab_high[sel])) +
  theme(legend.text = element_text(size=14))
ggsave(p2, file =paste0("ramkats_puiss_mean.pdf"), device =cairo_pdf, width= 16, height =4)

## 


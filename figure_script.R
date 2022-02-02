setwd("/mnt/9edae943-386a-4604-bb58-13212b0a3287/William/simulations_revisions/Figures/")
LowLinkage_PCst<-pcst_vals %>% filter(Linkage == "Low Linkage") %>% 
  mutate(Filter = fct_relevel(Filter, levels = "Out Combo", "Out Any", "Out Within","Out All","No Filter")) %>%
  ggplot(aes(x=V1,y=Filter))+  #theme_minimal() +
  geom_density_ridges(jittered_points = TRUE, alpha = 0.7,
                      point_size = 0.4, point_alpha = 0.5,rel_min_height = 0.000001,
                      quantile_lines=TRUE, 
                      quantile_fun=function(x,...)median(x),
                      vline_size = 1, vline_color = "red") +
  facet_grid(.~Scenario,scales="free") + theme_bw() + ylab("HWE Filter") + xlab("PCst")

ggsave(LowLinkage_PCst,filename = "LowLinkage_PCst.pdf",width = 11.6,height=3.86)

HighLinkage_PCst<-pcst_vals %>% filter(Linkage == "High Linkage") %>% 
  mutate(Filter = fct_relevel(Filter, levels = "Out Combo", "Out Any", "Out Within","Out All","No Filter")) %>%
  ggplot(aes(x=V1,y=Filter))+  #theme_minimal() +
  geom_density_ridges(jittered_points = TRUE, alpha = 0.7,
                      point_size = 0.4, point_alpha = 0.5,rel_min_height = 0.000001,
                      quantile_lines=TRUE, 
                      quantile_fun=function(x,...)median(x),
                      vline_size = 1, vline_color = "red") +
  facet_grid(.~Scenario,scales="free") + theme_bw() + ylab("HWE Filter") + xlab("PCst")

ggsave(HighLinkage_PCst,filename = "HighLinkage_PCst.pdf",width = 11.6,height=3.86)

LowLinkage_Fst <- fst_vals %>% filter(Linkage == "Low Linkage") %>%
  mutate(Filter = fct_relevel(Filter, levels = "Out Combo", "Out Any", "Out Within","Out All","No Filter")) %>%
  ggplot(aes(x=V1,y=Filter))+  #theme_minimal() +
  geom_density_ridges(jittered_points = TRUE, alpha = 0.7,
                      point_size = 0.4, point_alpha = 0.5,rel_min_height = 0.000001,
                      quantile_lines=TRUE, 
                      quantile_fun=function(x,...)median(x),
                      vline_size = 1, vline_color = "red") +
  facet_grid(.~Scenario,scales="free") + theme_bw() + ylab("HWE Filter") + xlab("Fst")

ggsave(LowLinkage_Fst, filename = "LowLinkage_Fst.pdf",width = 11.6,height=3.86)

HighLinkage_Fst <- fst_vals %>% filter(Linkage == "High Linkage") %>%
  mutate(Filter = fct_relevel(Filter, levels = "Out Combo", "Out Any", "Out Within","Out All","No Filter")) %>%
  ggplot(aes(x=V1,y=Filter))+  #theme_minimal() +
  geom_density_ridges(jittered_points = TRUE, alpha = 0.7,
                      point_size = 0.4, point_alpha = 0.5,rel_min_height = 0.000001,
                      quantile_lines=TRUE, 
                      quantile_fun=function(x,...)median(x),
                      vline_size = 1, vline_color = "red") +
  facet_grid(.~Scenario,scales="free") + theme_bw() + ylab("HWE Filter") + xlab("Fst")

ggsave(HighLinkage_Fst, filename = "HighLinkage_Fst.pdf",width = 11.6,height=3.86)

fis_regressions_plots_HighLinkage<-fis_regressions %>% filter(Linkage == "High Linkage") %>%
  mutate(Filter = fct_relevel(Filter, levels = "Out Combo", "Out Any", "Out Within","Out All","No Filter")) %>%
  ggplot(aes(x=Slope,y=Filter))+  #theme_minimal() +
  geom_density_ridges(jittered_points = TRUE, alpha = 0.7,
                      point_size = 0.4, point_alpha = 0.5,rel_min_height = 0.000001,
                      quantile_lines=TRUE, 
                      quantile_fun=function(x,...)median(x),
                      vline_size = 1, vline_color = "red") +
  facet_grid(.~Scenario,scales="fixed") + theme_bw() + ylab("HWE Filter") + xlab("Slope_OverallPop_Fis") 

ggsave(fis_regressions_plots_HighLinkage,filename = "HighLinkSlopes_FstVsOverallPopFis.pdf",width = 11.6,height=3.86)

fis_regressions_plots_LowLinkage<-fis_regressions %>% filter(Linkage == "Low Linkage") %>%
  mutate(Filter = fct_relevel(Filter, levels = "Out Combo", "Out Any", "Out Within","Out All","No Filter")) %>%
  ggplot(aes(x=Slope,y=Filter))+  #theme_minimal() +
  geom_density_ridges(jittered_points = TRUE, alpha = 0.7,
                      point_size = 0.4, point_alpha = 0.5,rel_min_height = 0.000001,
                      quantile_lines=TRUE, 
                      quantile_fun=function(x,...)median(x),
                      vline_size = 1, vline_color = "red") +
  facet_grid(.~Scenario,scales="fixed") + theme_bw() + ylab("HWE Filter") + xlab("Slope_OverallPop_Fis") 

ggsave(fis_regressions_plots_LowLinkage,filename = "LowLinkSlopes_FstVsOverallPopFis.pdf",width = 11.6,height=3.86)

HighLinkage_Struc<- k6_dat %>% filter(Linkage == "High Linkage") %>%
  mutate(filter = fct_relevel(filter, levels = "Out Combo", "Out Any", "Out Within","Out All","No Filter")) %>%
  ggplot(aes(x=k6_dat,y=filter))+  #theme_minimal() +
  geom_density_ridges(jittered_points = TRUE, alpha = 0.7,
                      point_size = 0.4, point_alpha = 0.5,rel_min_height = 0.000001,
                      quantile_lines=TRUE, 
                      quantile_fun=function(x,...)median(x),
                      vline_size = 1, vline_color = "red") +
  facet_grid(.~Scenario,scales="free") + theme_bw() + ylab("HWE Filter") + xlab("Average Nucleotide Distance")

ggsave(HighLinkage_Struc, filename = "HighLinkage_Struc.pdf",width = 11.6,height=3.86)

LowLinkage_Struc<- k6_dat %>% filter(Linkage == "Low Linkage") %>%
  mutate(filter = fct_relevel(filter, levels = "Out Combo", "Out Any", "Out Within","Out All","No Filter")) %>%
  ggplot(aes(x=k6_dat,y=filter))+  #theme_minimal() +
  geom_density_ridges(jittered_points = TRUE, alpha = 0.7,
                      point_size = 0.4, point_alpha = 0.5,rel_min_height = 0.000001,
                      quantile_lines=TRUE, 
                      quantile_fun=function(x,...)median(x),
                      vline_size = 1, vline_color = "red") +
  facet_grid(.~Scenario,scales="free") + theme_bw() + ylab("HWE Filter") + xlab("Average Nucleotide Distance")

ggsave(LowLinkage_Struc, filename = "LowLinkage_Struc.pdf",width = 11.6,height=3.86)



struc_real_plot<- struc_real %>% 
  mutate(filter = fct_relevel(filter, levels = "Out Combo", "Out Any", "Out Within","Out All","No Filter")) %>%
  mutate(Scen_F = fct_relevel(Scen_F, levels = "Seal","Zebra","Isopod")) %>%
  ggplot(aes(x=struc_mats,y=filter))+  #theme_minimal() +
  geom_density_ridges(jittered_points = TRUE, alpha = 0.7,
                      point_size = 0.4, point_alpha = 0.5,rel_min_height = 0.000001,
                      quantile_lines=TRUE, 
                      quantile_fun=function(x,...)median(x),
                      vline_size = 1, vline_color = "red") +
  facet_grid(.~Scen_F,scales="free") + theme_bw() + ylab("HWE Filter") + xlab("Average Nucleotide Distance")
ggsave(struc_real_plot, filename = "real_Struc.pdf",width = 11.6,height=3.86)



fst_real_plot<- real_fst %>% 
  mutate(Scenario = fct_relevel(Scenario, levels = "Out Combo", "Out Any", "Out Within","Out All","No Filter")) %>%
  mutate(DatType = fct_relevel(DatType, levels = "seal","zebra","isopod")) %>%
  ggplot(aes(x=V1,y=Scenario))+  #theme_minimal() +
  geom_density_ridges(jittered_points = TRUE, alpha = 0.7,
                      point_size = 0.4, point_alpha = 0.5,rel_min_height = 0.000001,
                      quantile_lines=TRUE, 
                      quantile_fun=function(x,...)median(x),
                      vline_size = 1, vline_color = "red") +
  facet_grid(.~DatType,scales="free") + theme_bw() + ylab("HWE Filter") + xlab("Fst")

ggsave(fst_real_plot, filename = "real_Fst.pdf",width = 11.6,height=3.86)



real_pcst_plot<-real_pcst %>% 
  mutate(Scenario = fct_relevel(Scenario, levels = "Out Combo", "Out Any", "Out Within","Out All","No Filter")) %>%
  mutate(DatType = fct_relevel(DatType, levels = "seal","zebra","isopod")) %>%
  ggplot(aes(x=V1,y=Scenario))+  #theme_minimal() +
  geom_density_ridges(jittered_points = TRUE, alpha = 0.7,
                      point_size = 0.4, point_alpha = 0.5,rel_min_height = 0.000001,
                      quantile_lines=TRUE, 
                      quantile_fun=function(x,...)median(x),
                      vline_size = 1, vline_color = "red") +
  facet_grid(.~DatType,scales="free") + theme_bw() + ylab("HWE Filter") + xlab("PCst") 

ggsave(real_pcst_plot,filename = "real_PCst.pdf",width = 11.6,height=3.86)


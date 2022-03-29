library(tidyverse)
library(rstan)
library(bayesplot)
library(scales) # to access break formatting functions
library(ggpmisc)
library(wesanderson)

load(here("data","Stan_Fit_combined_20210929_1338.RData"))

logtheta <- extract(Output$stanMod, par = "b_count")$b_count %>% colMeans()

 logtheta <- logtheta %>% 
  as.data.frame() %>% 
  mutate(station_id = Output$Stations$station_id) 
  
names(logtheta)[1:47] <- 1:47
pred <-  logtheta %>% 
   pivot_longer(-station_id, names_to = "sp_index_count") %>% 
   mutate(sp_index_count = as.numeric(sp_index_count)) 

Output$D_count %>% 
  dplyr::select(ID_count, Abund , station_id, sp_index_count) %>% 
  pivot_wider(., names_from="station_id", values_from ="Abund", values_fill=0) %>% 
  pivot_longer(cols=`1999_93.3_60`:`2013_86.7_50`, names_to="station_id", values_to ="Abund")-> Output$D_count_zero

 
obs <- Output$D_count_zero %>% 
  left_join(Output$species_mapping %>% dplyr::select(sp_index_main, sp_index_count)) %>% 
  left_join(Output$Stations) %>% 
  dplyr::select(station_id, station_idx, sp_index_count, Abund)

pred %>% 
  left_join(obs) %>% 
  drop_na() %>% 
  rename(predicted = value,
         observed = Abund) -> model_and_counts

model_and_counts %>% 
  ggplot(aes(x = observed, y = predicted))  + ylab("Predicted Counts") + xlab("Observed Morphological Counts") +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0,slope=1, linetype=2, color="tomato", size=2) +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),
                     limits=c(0, 10^5),
                     breaks=c(0,10^(1:5)))+
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10),
                     limits=c(0, 10^5),
                     breaks=c(0,10^(1:5))) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 14),
    panel.background = element_rect(fill = 'white', colour = 'white'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 14, angle = 0),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    strip.background = element_rect(fill = 'white'),
    strip.text.x = element_text(size = 15),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 15),
    axis.line = element_line(colour = 'black')
  ) + stat_cor(aes(label = ..r.label..),
               label.x = 1.5,
               label.y = 0.1, r.digits = 2,size = 7) -> model_vs_counts

ggsave(
  model_vs_counts,
  file = here::here("analysis", "figures", "S3_fig_model_vs_counts.png"),
  width = 12,
  height = 8
)


Output$D_mifish %>% 
  dplyr::select(station_id, nReads, sp_index_mifish) %>% 
  group_by(station_id, sp_index_mifish) %>% 
  summarize(nReads = sum(nReads)) %>% 
  pivot_wider(., names_from="station_id", values_from ="nReads", values_fill=0) %>% 
  pivot_longer(cols=`1996_80_60`:`2019_93.3_60`, names_to="station_id", values_to ="nReads") %>% 
  left_join(Output$species_mapping %>% dplyr::select(ID_main,sp_index_main, sp_index_count, sp_index_mifish)) %>% 
  right_join(Output$D_count_zero) -> reads_and_counts

reads_and_counts %>% 
  ggplot(aes(x = Abund, y = nReads))  + ylab("Observed Sequence Reads") +
  geom_point(alpha = 0.5) +
  xlab("Observed Morphological Counts") +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),
                     limits=c(0, 10^6),
                     breaks=c(0,10^(1:6)))+
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10),
                     limits=c(0, 10^4),
                     breaks=c(0,10^(1:4)))+
  theme_bw() +
  geom_abline(intercept = 0,slope=1, linetype=2, color="tomato", size=2) +
  theme(
    axis.text.x = element_text(size = 14),
    panel.background = element_rect(fill = 'white', colour = 'white'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 14, angle = 0),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    strip.background = element_rect(fill = 'white'),
    strip.text.x = element_text(size = 15),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 15),
    axis.line = element_line(colour = 'black')
  ) + stat_cor(aes(label = ..r.label..),
                   label.x = 1.5,
                   label.y = 0.1, r.digits = 2,size = 7) -> reads_vs_counts

ggsave(
  reads_vs_counts,
  file = here::here("analysis", "figures", "S2_fig_reads_vs_counts.png"),
  width = 12,
  height = 8
)

#posterior mean
po_1a <- Output$D_mifish %>% 
  ggplot(., aes(y = post_mean, x = nReads) ) +
  geom_point( data=Output$D_mifish, aes(y = post_mean, x = nReads), alpha = 0.2, color="black")  +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),
                     limits=c(0, 10^6),
                     breaks=c(0,10^(1:6)))+
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10),
                     limits=c(0, 10^6),
                     breaks=c(0,10^(1:6)))+
  geom_abline(intercept = 0,slope=1, linetype=2, color="tomato", size=2) +
  xlab("Observed Sequence Reads") +
  ylab("Predicted Sequence Reads") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 14),
    panel.background = element_rect(fill = 'white', colour = 'white'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 14, angle = 0),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    strip.background = element_rect(fill = 'white'),
    axis.line = element_line(colour = 'black')) +  stat_cor(aes(label = ..r.label..),
                                                             label.x = 1.5,
                                                             label.y = 0.1, r.digits = 2,size = 7)

po_1a

ggsave(
  po_1a,
  file = here::here("analysis", "figures", "S3_fig_pred_reads_vs_obs_reads.png"),
  width = 12,
  height = 8
)

# Supplemental Figure 5. Co-detection of eDNA and Microscopy Data


Output$D_mifish %>% 
  dplyr::select(station_id, nReads, sp_index_mifish) %>% 
  group_by(station_id, sp_index_mifish) %>% 
  summarize(nReads = sum(nReads)) %>% 
  pivot_wider(., names_from="station_id", values_from ="nReads", values_fill=0) %>% 
  pivot_longer(cols=`1996_80_60`:`2019_93.3_60`, names_to="station_id", values_to ="nReads") -> mifish_data_ready

Output$Stations -> station_holder 
Output$species_mapping -> species_mapping_holder
station_holder$fake <- 1
species_mapping_holder$fake <- 1
my_cross_join <- full_join(station_holder, species_mapping_holder, by = "fake") %>%
                     select(-fake)
                   
my_cross_join %>% 
  left_join(mifish_data_ready) %>% 
  left_join(Output$D_count_zero) %>% 
  dplyr::mutate(nReads = replace_na(nReads, 0),
                Abund = replace_na(Abund, 0))-> cross_all_jars


cross_all_jars %>% 
  mutate(., Morph_detect = ifelse(Abund > 0, 1,0)) %>% 
  mutate(., MiFish_detect = ifelse(nReads > 0, 1,0)) %>% 
  group_by(ID_main,station_id) %>% 
  mutate(type = case_when(MiFish_detect == 1 & Morph_detect == 1 ~"Both_detected",
                          MiFish_detect == 0 & Morph_detect == 0 ~"Both_absent",
                          MiFish_detect == 1 & Morph_detect == 0 ~"eDNA_only",
                          MiFish_detect == 0 & Morph_detect == 1 ~"Morph_only",
                          TRUE~"low_data")) %>% 
  mutate(., ID_main = str_replace(ID_main,"AANOTHER","sp.")) %>% 
  mutate(., ID_main = str_replace(ID_main,"nannochir","cold variant")) %>% 
  dplyr::select(Species=ID_main,station_id, MiFish_detect, Morph_detect, type) ->b_heat

wesanderson::wes_palette("Darjeeling1") -> col_dar
wesanderson::wes_palette("Darjeeling2") -> col_dar2

b_heat  %>% ggplot(., aes(x=Species, y=station_id, fill=type)) + 
  geom_tile() + theme(axis.text.x = element_text(angle = 90)) + scale_fill_manual(name = "Detection Type", labels =c("Both Absent","Both Detected","eDNA Only","Microscopy Only"),values = c("white",col_dar[2],col_dar[5],col_dar2[3])) +
  ylab("Year x Station")-> b_heat_plot

b_heat_plot


ggsave(
  b_heat_plot,
  file = here::here("analysis", "figures", "S5_fig_codetection_of_taxa.png"),
  width = 18,
  height = 18
)

b_heat %>% filter(., type=="Morph_only") %>% dplyr::select(ID_main=Species,station_id) %>% 
  left_join(cross_all_jars) %>% 
  arrange(desc(Abund))

cross_all_jars %>% 
  filter(., nReads <1) %>% 
  filter(., Abund >0) %>% arrange(desc(Abund)) %>% 
  dplyr::summarise(max(Abund),min(Abund),mean(Abund),sd(Abund))

b_heat %>% 
  group_by(Species) %>% 
  count(type) %>% 
  filter(.,type!="Both_absent") %>% 
  pivot_wider(names_from = type, values_from = n, values_fill=0) %>% 
  mutate(class = case_when(Both_detected>0 ~"Detected by Both at a site",
                            eDNA_only>0 & Morph_only >0 ~"Detected by Both but different sites",
                           eDNA_only>0 ~"eDNA Only",
                           Morph_only>0~ "Morph Only")) %>% 
  ungroup() %>% 
  group_by(class) %>% 
  count()


b_heat %>% 
  group_by(type) %>% 
  count() %>% 
  mutate( perc = n/4704)


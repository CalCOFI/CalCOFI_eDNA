library(tidyverse)
library(rstan)
library(bayesplot)
library(scales) # to access break formatting functions


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
  ) -> model_vs_counts

ggsave(
  model_vs_counts,
  file = here::here("analysis", "figures", "S1_fig_model_vs_counts.png"),
  width = 12,
  height = 8
)

Output$D_mifish %>% 
  dplyr::select(station_id, nReads, sp_index_mifish) %>% 
  group_by(station_id, sp_index_mifish) %>% 
  summarize(nReads = sum(nReads)) %>% 
  pivot_wider(., names_from="station_id", values_from ="nReads", values_fill=0) %>% 
  pivot_longer(cols=`1996_80_60`:`2019_93.3_60`, names_to="station_id", values_to ="nReads") %>% 
  left_join(Output$species_mapping %>% dplyr::select(sp_index_main, sp_index_count, sp_index_mifish)) %>% 
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
  ) -> reads_vs_counts

ggsave(
  reads_vs_counts,
  file = here::here("analysis", "figures", "S3_fig_reads_vs_counts.png"),
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
    axis.line = element_line(colour = 'black')) 

po_1a

ggsave(
  po_1a,
  file = here::here("analysis", "figures", "S2_fig_pred_reads_vs_obs_reads.png"),
  width = 12,
  height = 8
)

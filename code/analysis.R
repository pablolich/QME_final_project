library(tidyverse)
library(viridis)
library(RColorBrewer)

  setwd("~/Desktop/coursework/QME/QME_final_project/code/")
data = read.csv('../data/taylor_crossfeeding_results.csv')
data = read.csv('../data/cluster_results.csv')
data = read.csv('../data/cluster_results_try2.csv')


data2=read.csv('../data/cluster_results_try3.csv')
data2$sim = data2$sim + max(data$sim) + 1
data_comb = rbind(data, data2)

#average over simulations
data_w = data %>% 
  select(c(1:12, 14:17)) %>% 
  pivot_longer(c(distance_assembly, distance_optimal, distance_naive),
               names_to = "distance_method", values_to = "dist_value") %>% 
  group_by(leakage, n/N, distance_method) %>% 
  summarize(mean_dist = mean(dist_value), 
            sd_dist = sd(dist_value)) #%>% 
  #filter(mean_dist < 50)

ggplot(data_w , aes(x = `n/N`, y = log(1+mean_dist), 
                    group = interaction(distance_method, leakage), 
                    colour = as.factor(leakage)))+
  #geom_point() + 
  geom_line(aes(linetype = distance_method), size=1) +
  scale_linetype_manual(values=c('solid', 'dotted', 'longdash'))+
  theme_bw()+
  theme(legend.position = c(0.7, 0.7))+
  scale_color_brewer(palette = "YlOrRd")
  #scale_color_viridis_d("magma")
     #ggsave('../data/error_n_N.pdf', height = 4, width = 4.5)



library(tidyverse)
library(viridis)

setwd("~/Desktop/coursework/QME/QME_final_project/code/")
data = read.csv('../data/taylor_results.csv')

#average over simulations
data_w = data %>% 
  select(c(1:11, 13:15)) %>% 
  pivot_longer(c(distance_assembly, distance_optimal, distance_naive),
               names_to = "method", values_to = "dist_value") %>% 
  group_by(a, b, n/N, method) %>% 
  mutate(mean_dist = mean(dist_value))

ggplot(data_w , aes(x = n/N, y = mean_dist, 
                    group = method, 
                    colour = method))+
  geom_point() + 
  geom_line(aes(linetype = method)) +
  theme_bw()+
  theme(legend.position = c(0.7, 0.7))+
  ggsave('../data/error_n_N.pdf', height = 4, width = 4.5)



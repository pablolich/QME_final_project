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
data_w = data_comb %>% 
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
  scale_linetype_manual(values=c('solid', 'dotted', 'longdash'),
                        labels=c('Assembly', 'Null', 'Optimal'))+
  theme_bw()+
  theme(legend.position = c(0.8, 0.85),
        legend.box = "horizontal")+
  #scale_color_brewer(palette = "YlOrRd")
  scale_color_manual(values=c("#FF6B6B", "#FFD93D", "#6BCB77"))+
  labs(color="Leakage", linetype= "Method")
ggsave('../data/all_curves.pdf', width = 8.25, height = 6.50)

###############################################################################

#plot ROS of best subcommunity, summed accross simulations with the same N.
#each plot corresponds to a different value of leakage.

#Select portion of data to plot
data_ros = data_comb %>% select(c(X, sim, N, n, leakage, distance_assembly,
                                  ab_subcomm, 18:37))

dev.off()

preferences = data_ros %>% group_by(sim, N, n) %>% 
  filter(ab_subcomm != 0) %>% 
  slice_min(distance_assembly) %>% 
  pivot_longer(r1:r20, names_to = 'resource', values_to = 'use') %>% 
  group_by(sim, N, n, use, leakage) %>% 
  summarize(n_pref = mean(1/use)) %>% 
  filter(n_pref < Inf) %>%  
  group_by(N, n) %>% 
  mutate(n_mean = mean(n_pref), 
         n_var = sd(n_pref))

# The average number of preferences as a function of the size of the subcommunity 
# and the evolution of the mean vs varianvce of the distribution of generalists 
# and specialists as the number of species increases. This is always for the 
# subcommunity that best approximates the function of the original community. 
# The three panels are each for one value of leakage: 0, 0.5, 0.9. 
# The size of original community is N, and the size of the subcommunities is n. 
# I performed 100 simulations, and averaged over communities with the same N

ggplot(preferences, aes(x = n/N, y = n_mean, group = N,
                        color = N))+
  geom_point()+
  geom_line()+
  scale_color_viridis_c()+
  facet_wrap(~leakage, nrow = 3)+
  theme(aspect.ratio = 0.7)
ggsave('../data/average_n_preferences.pdf', width = 5.2, height = 10)

 preferences_filt <- preferences %>% filter(N >= 7)
ggplot(preferences_filt, aes(x = n_mean, y = n_var, fill = n/N))+
  geom_line(aes(group = interaction(N, sim), colour = N), 
            linetype = 'dashed', lwd = 0.5)+
  geom_point(shape = 21, lwd=2, color = 'grey')+
  facet_wrap(~leakage, nrow=3)+
  theme(aspect.ratio = 0.5)+
  scale_color_viridis()+ 
  scale_fill_viridis(option = 'mako')
ggsave('../data/mean_variance.pdf', height = 10, width = 7)

###############################################################################
  
#distribution of the best-subcommunity resource preferences, binned accross 
#simulations where the original community had the same diversity N
#different values of leakage are all lumped together here, otherwise we don't 
#have enough power.

pref_freq = preferences %>% 
  filter(N > 5) %>% 
  count(n, n_pref, name='freq') %>% 
  group_by(N, n) %>% 
  mutate(sum_freq = freq/sum(freq))
 ggplot(pref_freq)+
   geom_bar(aes(x=as.factor(n), y = sum_freq, fill = n_pref), 
            stat = 'identity')+
   scale_fill_viridis(option = 'turbo')+
   facet_wrap(~N, scales = 'free_x')+
   theme(aspect.ratio = 1)
ggsave('../data/preference_distribution_best.pdf', width = 11.9, height = 9)
 
 ##############################################################################
 
 #Plot worst community instead of best community
 preferences_worst = data_ros %>% group_by(sim, N, n) %>% 
   filter(ab_subcomm != 0) %>% 
   slice_max(distance_assembly) %>% 
   pivot_longer(r1:r20, names_to = 'resource', values_to = 'use') %>% 
   group_by(sim, N, n, use, leakage) %>% 
   summarize(n_pref = mean(1/use)) %>% 
   filter(n_pref < Inf) %>%  
   group_by(N, n) %>% 
   mutate(n_mean = mean(n_pref), 
          n_var = sd(n_pref))
 pref_freq_worst = preferences_worst %>% 
   filter(N > 5) %>% 
   count(n, n_pref, name='freq') %>% 
   group_by(N, n) %>% 
   mutate(sum_freq = freq/sum(freq))
 ggplot(pref_freq_worst)+
   geom_bar(aes(x=as.factor(n), y = sum_freq, fill = n_pref), 
            stat = 'identity')+
   scale_fill_viridis(option = 'turbo')+
   facet_wrap(~N, scales = 'free_x')+
   theme(aspect.ratio = 1)
ggsave('../data/preference_distribution_worst.pdf', width = 11.9, height = 9)
 
###############################################################################
 
#Plot communites of random species (taken from the pool of the assembled 
#subcomasd)
#1 For each N: Check the frequency of communities with richness N, f
#2 We know that for each replicate there is one subcommunity of each size. 
#Therefore, sample n species from each of the f communities of size N, for each
#value of n. Then store in a dataframe, and plot.

#get vector number samples per group
  
#build dataframe by sampling from each group with a different n

#split dataframe into groups
pref_split <- data_ros %>% 
  group_by(sim, N, n) %>% 
  filter(ab_subcomm != 0) %>% 
  group_split()
#get number of groups
n_groups = length(pref_split)
#Initialize dataframe to store each sample
pref_random = data_ros %>% filter(N == Inf)
#sample n elements from each group
for (i in seq(n_groups)){
  #sample n samples from each group
  n_samp = unique(pref_split[[i]]$n)
  sample_i = pref_split[[i]] %>% 
    slice_sample(n = n_samp)
  pref_random = rbind(pref_random, sample_i)
}
  
pref_freq_random <- pref_random %>% 
  pivot_longer(r1:r20, names_to = 'resource', values_to = 'use') %>% 
  group_by(sim, N, n, use, leakage) %>% 
  summarize(n_pref = mean(1/use)) %>% 
  filter(n_pref < Inf) %>%  
  group_by(N, n) %>% 
  mutate(n_mean = mean(n_pref), 
         n_var = sd(n_pref)) %>% 
  filter(N > 5) %>% 
  count(n, n_pref, name='freq') %>% 
  group_by(N, n) %>% 
  mutate(sum_freq = freq/sum(freq))
  

ggplot(pref_freq_random)+
  geom_bar(aes(x=as.factor(n), y = sum_freq, fill = n_pref), 
           stat = 'identity')+
  scale_fill_viridis(option = 'turbo')+
  facet_wrap(~N, scales = 'free_x')+
  theme(aspect.ratio = 1)
ggsave('../data/preference_distribution_random.pdf', width = 11.9, height = 9)


###############################################################################
#correlation between the abundances in the subcommunity and the abundances in
#the original community as n increases
abundances = data_comb %>% select(c(X,sim, N, n, leakage, distance_assembly,
                                  ab_subcomm, ab_original)) %>% 
  filter(ab_subcomm != 0) %>% 
  group_by(sim, N, n, leakage) %>% 
  slice(which(distance_assembly == min(distance_assembly))) %>% 
  mutate(corr = cor(ab_subcomm, ab_original)) %>% 
  ungroup(sim) %>% 
  mutate(av_corr = mean(corr))

ggplot(abundances, aes(x = n/N, y = av_corr))+
  geom_line(aes(group = N, color = N))+
  facet_wrap(~leakage, nrow=3)+
  scale_fill_viridis(option='mako')+
  scale_color_viridis()+
  theme(aspect.ratio = 0.7)
ggsave('../data/average_abundance_correlation.pdf', width = 5.2, height = 10)

###############################################################################

#are sub-communities self-contained within each other?
#Perform a wholistic analysis by looking at the nestedness of the network. 
#first, get the bipartite network
#2. Get the eigenvalues. if its nested I should see a biggest eigenvalue, 
#otherwise not. Read philip and alesina paper on the ghost of nestedness in 
#nature communications.

  


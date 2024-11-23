rm(list = ls())
# Generate the data
library(extraDistr)
set.seed(2024)
L.true = 5 # true number of groups in population 1
alpha0 = 10 
alpha1.true = rgamma(n = 1, shape = alpha0, rate = 1) # True alpha1 for generating data
alpha2.true = rgamma(n = 1, shape = alpha1.true, rate = 1) # True alpha2 for generating data
# True weights
beta1.true = as.numeric(rdirichlet(n = 1, alpha = rep(alpha1.true/L.true, L.true))) # True beta1
beta2.true = as.numeric(rdirichlet(n = 1, alpha = alpha2.true * beta1.true)) # True beta2
# Sample sizes
n1 = 100 # First population
n2 = 120 # Second population
z1.true = sample(1:L.true, size = n1, prob = beta1.true, replace = TRUE)
z2.true = sample(1:L.true, size = n2, prob = beta2.true, replace = TRUE)
# Draw data from Normal populations
mean.true = seq(from = -10, by = 7, len = L.true) # True mean
X1 = rnorm(n1, mean = mean.true[z1.true], sd = 1)  # Data from population 1
X2 = rnorm(n2, mean = mean.true[z2.true], sd = 1)  # Data from population 2 
library(tidyverse)
data.frame(sample = 1:length(X1), x = X1) %>% ggplot(aes(x = sample, y = x)) + geom_point() + labs(x = "Observation Index", y = "Y", title = "Population 1") + theme_classic() +  
  theme(
    # LABLES APPEARANCE
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    plot.title = element_text(hjust = 0.5, size=20, face= "bold", colour= "black" ),
    plot.subtitle = element_text(hjust = 0.5, size=16, face= "bold", colour= "black" ),
    axis.title.x = element_text(size=16, face="bold", colour = "black"),    
    axis.title.y = element_text(size=16, face="bold", colour = "black"),    
    axis.text.x = element_text(size=16, face="bold", colour = "black"), 
    axis.text.y = element_text(size=16, face="bold", colour = "black"),
    strip.text.x = element_text(size = 14, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 14, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.line.y = element_line(color="black", size = 0.3),
    panel.border = element_rect(colour = "black", fill=NA, size=0.3),
    legend.title=element_text(size=16),
    legend.text=element_text(size=14)
  ) 

data.frame(sample = 1:length(X2), x = X2) %>% ggplot(aes(x = sample, y = x)) + geom_point() + labs(x = "Observation Index", y = "Y", title = "Population 2") + theme_classic() +  
  theme(
    # LABLES APPEARANCE
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    plot.title = element_text(hjust = 0.5, size=20, face= "bold", colour= "black" ),
    plot.subtitle = element_text(hjust = 0.5, size=16, face= "bold", colour= "black" ),
    axis.title.x = element_text(size=16, face="bold", colour = "black"),    
    axis.title.y = element_text(size=16, face="bold", colour = "black"),    
    axis.text.x = element_text(size=16, face="bold", colour = "black"), 
    axis.text.y = element_text(size=16, face="bold", colour = "black"),
    strip.text.x = element_text(size = 14, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 14, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.line.y = element_line(color="black", size = 0.3),
    panel.border = element_rect(colour = "black", fill=NA, size=0.3),
    legend.title=element_text(size=16),
    legend.text=element_text(size=14)
  ) 

################################################################################
# Prior specifications
################################################################################
L.max = 10    # Truncation level of DP
alpha0 = 3    # Hyper-parameter of alpha1
lambda = 0.01 # Prior precision of Normal prior
xi = 0        # Prior mean of Normal prior

data_stan <- list(L = L.max,
                  n1 = n1,
                  n2 = n2,
                  
                  y1 = X1,
                  y2 = X2,
                  
                  alpha0 = alpha0,
                  xi = 0,
                  lambda = lambda
)


library(rstan)
time.start = Sys.time()
output <- stan(file = "/Users/arhitchakrabarti/Desktop/Coursework/Fall 2024/STAT 689/Codes/HMC_Sampling.stan", 
               data = data_stan, 
               # Below are optional arguments
               iter = 10000, 
               chains = 1, 
               warmup = 5000,
               thin = 5,
               cores = min(parallel::detectCores(), 4))
time.end = Sys.time()
time.end - time.start

round(unname(summary(output)$summary[2:11, 1]), 4)
post_est = summary(output)$summary[,1][-length(summary(output)$summary[,1])]
# Change the pars to only look at relevant outputs
# print(output, pars = "sigma")
# ShinyStan exploration
library(shinystan)
shinystan::launch_shinystan(output)
################################################################################
# Posterior probabilities from Stan output
post_prob1 <- matrix(summary(output, pars = "p1")$summary[, 1], ncol = L.max, byrow = TRUE)
post_prob2 <- matrix(summary(output, pars = "p2")$summary[, 1], ncol = L.max, byrow = TRUE)

cluster_assignment.post1 <- apply(post_prob1, MARGIN = 1, which.max)
cluster_assignment.post2 <- apply(post_prob2, MARGIN = 1, which.max)

ARI_pop1 = aricode::ARI(z1.true,
                        cluster_assignment.post1)
ARI_pop2 = aricode::ARI(z2.true,
                        cluster_assignment.post2)

# General clustering
cluster1 <- sort(unique(cluster_assignment.post1))
y1_cluster_index <- NULL
for(i in cluster1){
  y1_cluster_index <- rbind(y1_cluster_index, data.frame(observation = X1[cluster_assignment.post1 == i], cluster = i))
}

cluster2 <- sort(unique(cluster_assignment.post2))
y2_cluster_index <- NULL
for(i in cluster2){
  y2_cluster_index <- rbind(y2_cluster_index, data.frame(observation = X2[cluster_assignment.post2 == i], cluster = i))
}

y1_cluster_index$cluster <- as.factor(y1_cluster_index$cluster)
y2_cluster_index$cluster <- as.factor(y2_cluster_index$cluster)

library(tidyverse)
# Raw data clustering
plot1 <- data.frame(X1, cluster = as.factor(cluster_assignment.post1)) %>% ggplot(aes(x = 1:n1, y = X1, col = cluster)) + geom_point(size = 2) + labs(x = "Observation Index", y = "Y", title = "Population 1", subtitle = paste0("ARI = ", round(ARI_pop1, 5))) + theme_classic() +  
  theme(
    # LABLES APPEARANCE
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    plot.title = element_text(hjust = 0.5, size=20, face= "bold", colour= "black" ),
    plot.subtitle = element_text(hjust = 0.5, size=16, face= "bold", colour= "black" ),
    axis.title.x = element_text(size=16, face="bold", colour = "black"),    
    axis.title.y = element_text(size=16, face="bold", colour = "black"),    
    axis.text.x = element_text(size=16, face="bold", colour = "black"), 
    axis.text.y = element_text(size=16, face="bold", colour = "black"),
    strip.text.x = element_text(size = 14, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 14, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.line.y = element_line(color="black", size = 0.3),
    panel.border = element_rect(colour = "black", fill=NA, size=0.3),
    legend.title=element_text(size=16),
    legend.text=element_text(size=14)
  ) 

plot2 <- data.frame(X2, cluster = as.factor(cluster_assignment.post2)) %>% ggplot(aes(x = 1:n2, y = X2, col = cluster)) + geom_point(size = 2) + labs(x = "Observation Index", y = "Y", title = "Population 2", subtitle = paste0("ARI = ", round(ARI_pop2, 5))) + theme_classic() +  
  theme(
    # LABLES APPEARANCE
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    plot.title = element_text(hjust = 0.5, size=20, face= "bold", colour= "black" ),
    plot.subtitle = element_text(hjust = 0.5, size=16, face= "bold", colour= "black" ),
    axis.title.x = element_text(size=16, face="bold", colour = "black"),    
    axis.title.y = element_text(size=16, face="bold", colour = "black"),    
    axis.text.x = element_text(size=16, face="bold", colour = "black"), 
    axis.text.y = element_text(size=16, face="bold", colour = "black"),
    strip.text.x = element_text(size = 14, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 14, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.line.y = element_line(color="black", size = 0.3),
    panel.border = element_rect(colour = "black", fill=NA, size=0.3),
    legend.title=element_text(size=16),
    legend.text=element_text(size=14)
  ) 

gridExtra::grid.arrange(plot1, plot2, ncol = 2)

# Arranged data clustering
plot1_arranged <- y1_cluster_index %>% ggplot(aes(x = 1:n1, y = observation, col = cluster)) + geom_point(size = 2) + labs(x = "", y = "Observations", title = "Population 1", subtitle = paste0("ARI = ", round(ARI_pop1, 5))) + theme_classic() +  
  theme(
    # LABLES APPEARANCE
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    plot.title = element_text(hjust = 0.5, size=20, face= "bold", colour= "black" ),
    plot.subtitle = element_text(hjust = 0.5, size=16, face= "bold", colour= "black" ),
    axis.title.x = element_text(size=16, face="bold", colour = "black"),    
    axis.title.y = element_text(size=16, face="bold", colour = "black"),    
    axis.text.x = element_text(size=16, face="bold", colour = "black"), 
    axis.text.y = element_text(size=16, face="bold", colour = "black"),
    strip.text.x = element_text(size = 14, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 14, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.line.y = element_line(color="black", size = 0.3),
    panel.border = element_rect(colour = "black", fill=NA, size=0.3),
    legend.title=element_text(size=16),
    legend.text=element_text(size=14)
  ) 

plot2_arranged <- y2_cluster_index %>% ggplot(aes(x = 1:n2, y = observation, col = cluster)) + geom_point(size = 2) + labs(x = "", y = "Observations", title = "Population 2", subtitle = paste0("ARI = ", round(ARI_pop2, 5))) + theme_classic() +  
  theme(
    # LABLES APPEARANCE
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    plot.title = element_text(hjust = 0.5, size=20, face= "bold", colour= "black" ),
    plot.subtitle = element_text(hjust = 0.5, size=16, face= "bold", colour= "black" ),
    axis.title.x = element_text(size=16, face="bold", colour = "black"),    
    axis.title.y = element_text(size=16, face="bold", colour = "black"),    
    axis.text.x = element_text(size=16, face="bold", colour = "black"), 
    axis.text.y = element_text(size=16, face="bold", colour = "black"),
    strip.text.x = element_text(size = 14, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 14, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.line.y = element_line(color="black", size = 0.3),
    panel.border = element_rect(colour = "black", fill=NA, size=0.3),
    legend.title=element_text(size=16),
    legend.text=element_text(size=14)
  ) 

gridExtra::grid.arrange(plot1_arranged, plot2_arranged, ncol = 2)







cluster_STAN <- data.frame(Index = c(1:n1, 1:n2),
                           Y = c(X1, X2),
                           Cluster = factor(c(cluster_assignment.post1,
                                              cluster_assignment.post2)),
                           Population = factor(c(rep("Population 1", length(X1)),
                                                 rep("Population 2", length(X2)))))

cluster_STAN$Population_ARI <- factor(c(rep(paste0("Population 1\nARI = ", round(ARI_pop1, 5)), length(X1)),
                                        rep(paste0("Population 2\nARI = ", round(ARI_pop2, 5)), length(X2))
))


myvalues_color =  c("1" = "#F8766D",
                    "2" = "#00BA38",
                    "3" = "#619CFF", 
                    "4" = "blueviolet",
                    "5" = "cyan4",
                    "6" = "#E6AB02",
                    "7" = "#E36EF6",
                    "8" = "bisque4",
                    "9" = "coral4",
                    "10" = "darkslateblue",
                    "11" = "lightseagreen",
                    "12" = "#E69F00", 
                    "13" = "#AA3377",
                    "14" = "sienna3",
                    "15" = "hotpink",
                    "16" = "sienna4",
                    "17" = "hotpink3",
                    "18" = "sienna1",
                    "19" = "dodgerblue4",
                    "20" = "bisque2")

library(patchwork)
library(tidyverse)
library(latex2exp)

x.limit.lower = 1
x.limit.upper = max(n1, n2)

y.limit.lower = min(c(X1, X2))
y.limit.upper = max(c(X1, X2))

plot.clustering_STAN <- cluster_STAN %>%
  ggplot(aes(x = Index, y = Y, col = Cluster)) + geom_point(size = 2) + labs(x = "Observation Index", y = TeX("$Y$")) + facet_wrap(~Population_ARI, ncol = 2) + 
  scale_color_manual(values = myvalues_color) + xlim(x.limit.lower, x.limit.upper) + ylim(y.limit.lower, y.limit.upper) + theme_classic() +  
  theme(
    # LABLES APPEARANCE
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    plot.title = element_text(hjust = 0.5, size=20, face= "bold", colour= "black" ),
    plot.subtitle = element_text(hjust = 0.5, size=16, face= "bold", colour= "black" ),
    axis.title.x = element_text(size=16, face="bold", colour = "black"),    
    axis.title.y = element_text(size=16, face="bold", colour = "black"),    
    axis.text.x = element_text(size=16, face="bold", colour = "black"), 
    axis.text.y = element_text(size=16, face="bold", colour = "black"),
    strip.text.x = element_text(size = 14, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 14, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.line.y = element_line(color="black", size = 0.3),
    panel.border = element_rect(colour = "black", fill=NA, size=0.3),
    legend.title=element_text(size=16),
    legend.text=element_text(size=14)
  ) 

plot.clustering_STAN

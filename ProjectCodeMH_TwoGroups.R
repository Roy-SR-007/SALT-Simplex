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
L.max = 10 # Truncation level of DP
alpha0 = 3
alpha1.start = rgamma(n = 1, shape = alpha0) # Prior for Dirichlet distribution for group 1
alpha2.start = rgamma(n = 1, shape = alpha1.start) # Prior for Dirichlet distribution for group 2

xi = 0 # Prior mean of Normal
lambda = 0.01 # Prior precision of Normal

mu.start = rnorm(n = L.max, mean = xi, sd = 1/sqrt(lambda))
beta1.start = as.numeric(rdirichlet(n = 1, alpha = rep(alpha1.start/L.max, L.max))) # This becomes the starting value of beta1 in the Gibbs iteration
# To avoid numerical issues I add a small quantity to each coordinate and re-normalize
beta1.start = beta1.start + 1e-8 # Add 10^-8
beta1.start = beta1.start/sum(beta1.start) 
# Draw starting values of latent indicators for group 1
z1.start = extraDistr::rcat(n = n1, prob = beta1.start)

## Starting values of second population
beta2.start = as.numeric(rdirichlet(n = 1, alpha = alpha2.start * beta1.start))
# To avoid numerical issues I add a small quantity to each coordinate and re-normalize
beta2.start = beta2.start + 1e-8 # Add 10^-8
beta2.start = beta2.start/sum(beta2.start) # Re-normalize
# Draw starting values of latent indicators for group 2
z2.start = extraDistr::rcat(n = n2, prob = beta2.start)
#####################################################################################################
## Initialize storage
#####################################################################################################
z1_samples <- list(); z2_samples <- list()  # Latent indicators

beta1_samples <- list() ; beta2_samples <- list() # Weights
n1_samples <- list() ; n2_samples <- list() # To save cluster specific sample sizes
phi_samples <- list() # Storing the means of the normal distribution
alpha1_samples <- list(); alpha2_samples <- list() # To store the concentration parameters of the GDP
prob1_samples <- list(); prob2_samples <- list() # To store the probabilities of the two groups
# Initialize all the values at the starting values
beta1_samples[[1]] <- beta1.start; beta2_samples[[1]] <- beta2.start
phi_samples[[1]] <- mu.start 
z1_samples[[1]] <- z1.start; z2_samples[[1]] <- z2.start 
alpha1_samples[[1]] <- alpha1.start; alpha2_samples[[1]] <- alpha2.start 

accept_beta1 = 0
#####################################################################################################
# Functions
#####################################################################################################
# Write a function which gives mean is given vector is non-null and returns zero if the vector is null
mean_my <- function(x){
  if(length(x) == 0){
    return(0)
  }else{
    return(mean(x))
  }
}

sample_z <- function(data, beta, phi, L, dense = TRUE){
  n = length(data)
  # Define a matrix to calculate the probabilities
  probability = matrix(0, nrow = n, ncol = L)
  z = 0 # Define a vector of indicators
  for(i in 1:n){
    for(h in 1:L){
      probability[i, h] = log(beta[h]) + dnorm(data[i], mean = phi[h], log = TRUE)
    }
    # Use the log sum exp trick 
    probability[i, ] = exp(probability[i, ] - max(probability[i, ])) 
    # Return the probability
    probability[i, ] = probability[i, ]/sum(probability[i, ])
    z[i] = extraDistr::rcat(n = 1, prob = probability[i, ])
  }
  if(dense == TRUE){
    z = dplyr::dense_rank(z)
  }
  
  return(list(z = z, probability = probability))
}

# Define the target function for [beta1| - ] in log scale with Metropolis ratio
# @pars argument should be a named list with alpha1, alpha2, beta2
# @ycand corresponds to a proposed value of beta1
# @n1 corresponds to n1_h 's at any stage
log_beta1 <- function(ycand, n1, pars) {
  L = length(n1)
  alpha1 = as.numeric(pars$alpha1); alpha2 = as.numeric(pars$alpha2); beta2 = pars$beta2
  
  out = sum((n1 + alpha1/L - 1)*log(ycand) + (alpha2*ycand - 1)*log(beta2) - lgamma(alpha2 * ycand))
  return(out)
}  

# @pars argument should be a named list with alpha0, alpha2, beta1
# @ycurrent corresponds to the previous value of alpha1
sample_alpha1 <- function(ycurrent, pars){
  # Extract the parameters from the pars dataframe
  beta1 = as.numeric(pars$beta1)
  alpha0 = pars$alpha0
  alpha2 = pars$alpha2
  L = length(beta1) # Length of the simplex
  
  ycand = rgamma(n = 1, shape = alpha0) # Propose a candidate
  
  # Numerator of MH ratio in log-scale
  num = sum(ycand/L * log(beta1)) - L*lgamma(ycand/L) - ycand + (alpha0 - 1)*log(ycand) + ycand * log(alpha2) + dgamma(ycurrent, shape = alpha0, log = TRUE)
  
  # Denominator of MH ratio in log-scale
  den = sum(ycurrent/L * log(beta1)) - L*lgamma(ycurrent/L) - ycurrent + (alpha0 - 1)*log(ycurrent) + ycurrent * log(alpha2) + dgamma(ycand, shape = alpha0, log = TRUE)
  
  
  MH = num - den # MH ratio which become difference in log-scale
  if(log(runif(n = 1)) < MH){
    out = ycand # Accept proposed value if log(U) < MH
    accept = 1
  }else{
    out = ycurrent # Reject proposed value if log(U) < MH
    accept = 0
  }
  return(list(out = out, accept = accept))
}

# @pars argument should be a named list with alpha1, beta1, beta2
# @ycurrent corresponds to the previous value of alpha4
sample_alpha2 <- function(ycurrent, pars){
  # Extract the parameters from the pars dataframe
  beta1 = as.numeric(pars$beta1); beta2 = as.numeric(pars$beta2)
  alpha1 = pars$alpha1
  L = length(beta1) # Length of the simplex
  
  ycand = rgamma(n = 1, shape = alpha1) # Propose a candidate
  # Numerator of MH ratio in log-scale
  num = lgamma(ycand) + sum((ycand * beta1 * log(beta2)) - lgamma(ycand * beta1)) - ycand + (alpha1 - 1)*log(ycand) + dgamma(ycurrent, shape = alpha1, log = TRUE)
  # Denominator of MH ratio in log-scale
  den = lgamma(ycurrent) + sum((ycurrent * beta1 * log(beta2)) - lgamma(ycurrent * beta1)) - ycurrent + (alpha1 - 1)*log(ycurrent) + dgamma(ycand, shape = alpha1, log = TRUE)
  
  MH = num - den # MH ratio which become difference in log-scale
  if(log(runif(n = 1)) < MH){ 
    out = ycand # Accept proposed value if log(U) < MH
    accept = 1
  }else{
    out = ycurrent # Reject proposed value if log(U) < MH
    accept = 0
  }
  return(list(out = out, accept = accept))
}

indicator <- function(x, l){
  if(x == l){
    return(1)
  }else{
    return(0)
  }
}
################################################################################
# Running the Sampler
################################################################################
L.obs = L.max
gibbs_iterations = 50000
prior_prec = lambda # Prior mean of Normal distribution of G0
prior_mean = xi # Prior precision of Normal distribution of G0
rho = 0.99 # Tuning parameter of Dirichlet proposal distribution


for(b in 2:gibbs_iterations){
  # b = 2
  # Printing the iterations
  if(b == 2){
    cat(paste0("Iteration: ", (b-1), "\n"))
  }
  if(b %% floor((10/100)*(gibbs_iterations + 1)) == 0) {
    cat(paste0("Iteration: ", b, "\n"))
  }
  
  # Sample the latent indicators for population 1
  sample1 = sample_z(data = X1, beta = beta1_samples[[b-1]], phi = phi_samples[[b-1]], L = L.obs, dense = FALSE)
  z1_samples[[b]] = sample1$z# Store the samples values of z1
  # Sample the latent indicators for population 2
  sample2 = sample_z(data = X2, beta = beta2_samples[[b-1]], phi = phi_samples[[b-1]], L = L.obs, dense = FALSE)
  z2_samples[[b]] = sample2$z # Store the samples values of z2
  
  n1_h = 0; n2_h = 0 # Initialize n_h
  x1_h = 0; x2_h = 0 # Initialize mean_h
  phi = 0
  # L denotes the number of draws from the normal distribution of Full conditional
  # Note: n.draw should be same as the number of clusters as they correspond to an update for each cluster
  
  for(h in 1:L.obs){
    n1_h[h] = sum(z1_samples[[b]] == h) # Cluster size for population 1
    n2_h[h] = sum(z2_samples[[b]] == h) # Cluster size for population 2
    x1_h[h] = mean_my(X1[z1_samples[[b]] == h]) # Cluster mean for population 1
    x2_h[h] = mean_my(X2[z2_samples[[b]] == h]) # Cluster mean for population 2
    # Draw means
    phi[h] = rnorm(n = 1, 
                   mean = (n1_h[h] * x1_h[h] + n2_h[h] * x2_h[h] + prior_prec*prior_mean)/(prior_prec + n1_h[h] + n2_h[h]), 
                   sd = 1/sqrt((prior_prec + n1_h[h] + n1_h[h]))
    )
  }
  # Store the samples of phi, n1 and n2
  phi_samples[[b]] = phi; n1_samples[[b]] = n1_h; n2_samples[[b]] = n2_h
  
  # Run MH for Beta1
  # beta1.proposed = as.numeric(rdirichlet(n = 1, alpha = rep(rho * alpha1_samples[[b-1]]/L.max, L.max)))
  # beta1.proposed = as.numeric(rdirichlet(n = 1, alpha = rep(rho * 1/L.max, L.max)))
  beta1.proposed = as.numeric(rdirichlet(n = 1, alpha = rho * (alpha1_samples[[b-1]]/L.max + n1_samples[[b]])))
  beta1.proposed = beta1.proposed + 1e-8
  beta1.proposed = beta1.proposed/sum(beta1.proposed)
  
  log_num = log_beta1(ycand = beta1.proposed, n1 = n1_samples[[b]], pars = list(beta2 = beta2_samples[[b-1]],
                                                                            alpha1 = alpha1_samples[[b-1]],
                                                                            alpha2 = alpha2_samples[[b-1]])) +
    # extraDistr::ddirichlet(x = beta1_samples[[b-1]], alpha = rep(rho * alpha1_samples[[b-1]]/L.max, L.max), log = TRUE)
    extraDistr::ddirichlet(x = beta1_samples[[b-1]], alpha = rho * (alpha1_samples[[b-1]]/L.max + n1_samples[[b]]), log = TRUE)
    # extraDistr::ddirichlet(x = beta1_samples[[b-1]], alpha = rep(rho * 1/L.max, L.max), log = TRUE)
  
  log_den = log_beta1(ycand = beta1_samples[[b-1]], n1 = n1_samples[[b]], pars = list(beta2 = beta2_samples[[b-1]],
                                                                            alpha1 = alpha1_samples[[b-1]],
                                                                            alpha2 = alpha2_samples[[b-1]])) +
    # extraDistr::ddirichlet(x = beta1.proposed, alpha = rep(rho * alpha1_samples[[b-1]]/L.max, L.max), log = TRUE)
    extraDistr::ddirichlet(x = beta1.proposed, alpha = rho * (alpha1_samples[[b-1]]/L.max + n1_samples[[b]]), log = TRUE)
    # extraDistr::ddirichlet(x = beta1.proposed, alpha = rep(rho * 1/L.max, L.max), log = TRUE)
  
  log_alpha = min(0, (log_num - log_den))
  u = runif(n = 1)
  
  beta1_samples[[b]] = (log(u) < log_alpha)*beta1.proposed + (log(u) >= log_alpha)*beta1_samples[[b-1]]
  accept_beta1[b] = (log(u) < log_alpha)*1 + (log(u) >= log_alpha)*0

  # Draw samples from the beta2
  new_concentration = alpha2_samples[[b-1]] * beta1_samples[[b]]
  # Now we only need to consider components 1, 2,...,L.obs 
  beta2_samples[[b]] = as.numeric(rdirichlet(n = 1, alpha = n2_samples[[b]] + new_concentration))
  # To avoid numerical issues I add a small quantity to each coordinate and re-normalize
  beta2_samples[[b]] = beta2_samples[[b]] + 1e-8 
  beta2_samples[[b]] = beta2_samples[[b]]/sum(beta2_samples[[b]]) # Re-normalize
  
  alpha1_samples[[b]] <- sample_alpha1(ycurrent = alpha1_samples[[b-1]], 
                                       pars = list(alpha0 = alpha0, 
                                                   alpha2 = alpha2_samples[[b-1]], 
                                                   beta1 = beta1_samples[[b]]))$out
  
  alpha2_samples[[b]] <- sample_alpha2(ycurrent = alpha2_samples[[b-1]], 
                                       pars = list(alpha1 = alpha1_samples[[b]],
                                                   beta1 = beta1_samples[[b]],
                                                   beta2 = beta2_samples[[b]]))$out
  
}

gibbs_burn <- 25000
thin = 25
thin_samples <- seq(from = (gibbs_burn + 1), to = gibbs_iterations, by = thin)
mean(accept_beta1[thin_samples])

ll <- 0
for(b in 2:gibbs_iterations){
  # Printing the iterations
  if(b == 2){
    cat(paste0("Iteration: ", (b-1), "\n"))
  }
  if(b %% floor((10/100)*(gibbs_iterations + 1)) == 0) {
    cat(paste0("Iteration: ", b, "\n"))
  }
  
  log_like1 = 0; log_like2 = 0
  for(i in 1:n1){
    mixture1 = 0
    for(l in 1:L.max){
      mixture1[l] = indicator(z1_samples[[b]][i], l) * (log(beta1_samples[[b]][l]) + dnorm(X1[i], mean =  phi_samples[[b]][[l]], sd = 1, log = TRUE))
    }
    log_like1[i] = sum(mixture1)
  }
  
  for(i in 1:n2){
    mixture2 = 0
    for(l in 1:L.max){
      mixture2[l] = indicator(z2_samples[[b]][i], l) * (log(beta2_samples[[b]][l]) + dnorm(X2[i], mean =  phi_samples[[b]][[l]], sd = 1, log = TRUE))
    }
    log_like2[i] = sum(mixture2)
  }
  ll[b] = sum(c(log_like1, log_like2))
  
}

gibbs_burn <- 25000
thin = 25
thin_samples <- seq(from = (gibbs_burn + 1), to = gibbs_iterations, by = thin)
#####################################################################################
## Convergence and Mixing
#####################################################################################
ll_post_burn <- ll[thin_samples]

mixture1.true = 0; log_like1.true = 0
for(i in 1:n1){
  for(h in 1:L.true){
    mixture1.true[h] = indicator(z1.true[i], h) * (log(beta1.true[h]) + dnorm(X1[i], mean = mean.true[h], sd = 1, log = TRUE))
  }
  log_like1.true[i] = sum(mixture1.true)
}
mixture2.true = 0; log_like2.true = 0
for(i in 1:n2){
  for(h in 1:L.true){
    mixture2.true[h] = indicator(z2.true[i], h) * (log(beta2.true[h]) + dnorm(X2[i], mean = mean.true[h], sd = 1, log = TRUE))
  }
  log_like2.true[i] = sum(mixture2.true)
}

ll.true = sum(c(log_like1.true, log_like2.true))

ll_plot <- data.frame(x = 1:length(ll_post_burn), y = ll_post_burn) %>% ggplot(aes(x = x, y = y)) + geom_line() +
  labs(x = "Iteration post burn-in", y = "") +
  theme_classic() +  
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
    panel.border = element_rect(colour = "black", fill=NA, size=0.3)
  ) + geom_hline(yintercept = ll.true,
                 color = "red")

library(LaplacesDemon)
ll_plot <- ll_plot + labs(subtitle = paste0("Geweke Diagnostic Statistic = ", round(Geweke.Diagnostic(matrix(ll_post_burn, ncol = 1)), 4)))

library(forecast)
ACF_plot <-  ggAcf(x = ll_post_burn, lag.max = 40) + labs(title = "", y = "") +
  theme_classic() +  
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
    panel.border = element_rect(colour = "black", fill=NA, size=0.3)
  )

gridExtra::grid.arrange(ll_plot, ACF_plot, ncol = 2)

library(latex2exp)
alpha1_plot <- data.frame(alpha = unlist(alpha1_samples[thin_samples]),
                          Iteration = 1:length(thin_samples)) %>% ggplot(aes(x = Iteration, y = alpha)) + geom_line() +
  labs(
    x = "Iteration post burn-in", y = TeX("$\\alpha_1$")) +
  theme_classic() +
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


alpha2_plot <- data.frame(alpha = unlist(alpha2_samples[thin_samples]),
                          Iteration = 1:length(thin_samples)) %>% ggplot(aes(x = Iteration, y = alpha)) + geom_line() +
  labs(
    x = "Iteration post burn-in", y = TeX("$\\alpha_2$")) +
  theme_classic() +
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

gridExtra::grid.arrange(alpha1_plot, alpha2_plot, ncol = 2)

beta1_samples_matrix = matrix(unlist(beta1_samples[thin_samples]), nrow = length(thin_samples), ncol = L.max, byrow = TRUE)

plot_beta1 <- list()
for(i in 1:L.max){
  plot_beta1[[i]] <- data.frame(beta = beta1_samples_matrix[, i],
                                Iteration = 1:length(thin_samples)) %>% ggplot(aes(x = Iteration, y = beta)) + geom_line() +
    labs(
      x = "Iteration post burn-in", y = TeX(paste0("$\\beta_{1", i, "}$"))) +
    theme_classic() +
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
}

library(gridExtra)
grid.arrange(grobs = plot_beta1, ncol = 5)
#####################################################################################
## Clustering by mimimizing Variation of Information
#####################################################################################
library(mcclust.ext)

latent_indicators_combined <- list()
for(b in 2:gibbs_iterations){
  latent_indicators_combined[[b]] = c(z1_samples[[b]],
                                      z2_samples[[b]])
}

posterior_samples_latent_indicators_combined <- matrix(0, nrow = length(thin_samples), ncol = length(latent_indicators_combined[thin_samples][[1]]))

for(i in 1:length(thin_samples)){
  posterior_samples_latent_indicators_combined[i, ] = latent_indicators_combined[thin_samples][[i]]
}

sim_mat = comp.psm(posterior_samples_latent_indicators_combined)
par(mfrow = c(1, 1))
plotpsm(sim_mat)

################################################################################
# CLUSTERING BY MINIMIZING VI
################################################################################
best_cluster = minVI(sim_mat, posterior_samples_latent_indicators_combined, method = "all")

summary(best_cluster)

z.estimated_combined = best_cluster$cl[1, ]

cluster_VI <- data.frame(Index = c(1:n1, 1:n2),
                         Y = c(X1, X2),
                         Cluster = factor(z.estimated_combined),
                         Population = factor(c(rep("Population 1", length(X1)),
                                               rep("Population 2", length(X2)))))

ARI_pop1 = aricode::ARI(z1.true,
                        cluster_VI %>% filter(Population == "Population 1") %>% pull(Cluster))
ARI_pop2 = aricode::ARI(z2.true,
                        cluster_VI %>% filter(Population == "Population 2") %>% pull(Cluster))

cluster_VI$Population_ARI <- factor(c(rep(paste0("Population 1\nARI = ", round(ARI_pop1, 5)), length(X1)),
                                      rep(paste0("Population 2\nARI = ", round(ARI_pop2, 5)), length(X2))
)
)
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

plot.clustering <- cluster_VI %>%
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

plot.clustering

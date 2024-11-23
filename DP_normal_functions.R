################################################################################
# Functions
################################################################################
sample_mu <- function(data, z, L, prior_mean, prior_prec = 1){
  n_h = 0 # Initialize n_h
  x_h = 0 # Initialize mean_h
  n_h = 0; x_h = 0; phi = 0
  # L denotes the number of draws from the normal distribution of Full conditional
  # Note: n.draw should be same as the number of clusters as they correspond to an update for each cluster
  
  for(h in 1:L){
    n_h[h] = sum(z == h) # Cluster size
    x_h[h] = mean(data[z == h]) # Cluster mean
    phi[h] = rnorm(n = 1, mean = (n_h[h] * x_h[h] + prior_prec*prior_mean)/(prior_prec + n_h[h]), 
                   sd = 1/sqrt((prior_prec + n_h[h]))
    )
  }
  return(list(phi = phi, n_h = n_h))
}

sample_z <- function(data, beta, phi, L){
  n = length(data)
  # Define a matrix to calculate the probabilities
  probability = matrix(0, nrow = n, ncol = L)
  z = 0 # Define a vector of indicators
  for(i in 1:n){
    for(h in 1:L){
      # probability[i, h] = beta[h] * dnorm(data[i], mean = phi[h])
      probability[i, h] = log(beta[h]) + dnorm(data[i], mean = phi[h], log = TRUE)
    }
    # Use the log sum exp trick 
    probability[i, ] = exp(probability[i, ] - max(probability[i, ])) 
    # Return the probability
    probability[i, ] = probability[i, ]/sum(probability[i, ])
    z[i] = extraDistr::rcat(n = 1, prob = probability[i, ])
  }
  
  z = dplyr::dense_rank(z)
  return(list(z = z, probability = probability))
}

gibbs_dp <- function(init, data, num_iterations, burn, thin = 1, L.max, prior, samples = FALSE){
  n = length(data)
  # Define the lists to store the samples
  beta_samples <- list() # list to store samples of weights beta
  phi_samples <- list() # list to store Normal means or phi samples 
  z_samples <- list() # list to store Normal latent indicator samples 
  n_samples <- list() # list to store latent cluster sizes
  probability_samples <- list() # list to store matrix of probabilities
  
  # Initialize the starting values of the different parameters
  beta_samples[[1]] <- init$beta.start
  phi_samples[[1]] <- init$mu.start 
  z_samples[[1]] <- init$z.start 
  
  # Initialize L = L.max
  L.obs = L.max
  
  for(i in 2:num_iterations){
    # Printing the iterations
    if(i == 2){
      cat(paste0("Iteration: ", (i-1), "\n"))
    }
    if(i %% floor((10/100)*(num_iterations + 1)) == 0) {
      cat(paste0("Iteration: ", i, "\n"))
    }
    # Sample the latent cluster indicators and matrix of probabilities
    sample = sample_z(data = data, beta = beta_samples[[i-1]], phi = phi_samples[[i-1]], L = L.obs)
    
    z_samples[[i]] = sample$z # Take out the latent cluster indicators
    probability_samples[[i]] = sample$probability # Take out the matrix of probabilities
    
    L.obs = max(z_samples[[i]]) # Update L based on new cluster labels. The maximum truncation level of the DP is changed to the maximum value of the cluster indicators
    # Sample the means of the normal dsitribution
    phi_sampling = sample_mu(data = data, z = z_samples[[i]], L = L.obs, prior_mean = prior$prior_mean, prior_prec = prior$prior_prec)
    
    phi_samples[[i]] = phi_sampling$phi # Means
    n_samples[[i]] = phi_sampling$n_h # Cluster sizes
    
    # Sample the weights
    beta_samples[[i]] = as.numeric(rdirichlet(n = 1, alpha =   n_samples[[i]] + (prior$alpha)/L.obs))
    
  }# End of Gibbs update
  # Define posterior matrix to store the posterior samples post burn-in and after thinning
  phi_samp_post = NULL; beta_samp_post = NULL
  
  thin_samples <- seq(from = (burn + 1), to = num_iterations, by = thin)
  
  for(b in thin_samples){
    phi_samp_post = rbind(phi_samp_post, phi_samples[[b]])
    beta_samp_post = rbind(beta_samp_post, beta_samples[[b]])
  }
  
  if(samples == TRUE){
    return(list(phi_samp_post = phi_samp_post,
                beta_samp_post = beta_samp_post,
                probability_samples_all = probability_samples))
  }else{
    return(list(phi_post_mean = colMeans(phi_samp_post),
                beta_post_mean = colMeans(beta_samp_post)))
  }
}
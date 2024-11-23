// The input data is a vector 'y1' of length 'n1' for group 1 
// The input data is a vector 'y2' of length 'n2' for group 2
// Maximum number of mixture components 'L'
data {
  int<lower=1> L;           // number of mixture components
  int<lower=1> n1;           // number of data points in group 1
  int<lower=1> n2;           // number of data points in group 2
  vector[n1] y1;             // observations in group 1
  vector[n2] y2;             // observations in group 2

  real xi;
  real<lower=0> lambda;
  real<lower=0> alpha0;
}

parameters {
  real<lower=0> alpha1;      // concentration parameter of population 1
  simplex[L] beta1;          // mixing proportions 1
  ordered[L] mu;             // locations of mixture components 
  
  real<lower=0> alpha2;      // concentration parameter of population 2
  simplex[L] beta2;          // mixing proportions 2

}


// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  alpha1 ~ gamma(alpha0, 1); // Gamma prior with known alpha0
  beta1 ~ dirichlet(rep_vector((alpha1/L), L));
  alpha2 ~ gamma(alpha1, 1); // Gamma prior conditional on alpha1
  
  beta2 ~ dirichlet(alpha2 * beta1);

  mu ~ normal(rep_vector(xi, L), (1/lambda));

  for (n in 1:n1) {
    vector[L] lps1 = log(beta1);
    for (l in 1:L) {
      lps1[l] += normal_lpdf(y1[n] | mu[l], 1);
    }
    target += log_sum_exp(lps1);
  }
  
  for (n in 1:n2) {
    vector[L] lps2 = log(beta2);
    for (l in 1:L) {
      lps2[l] += normal_lpdf(y2[n] | mu[l], 1);
    }
    target += log_sum_exp(lps2);
  }
}

generated quantities {
  matrix[n1, L] p1;
  matrix[n2, L] p2;
  
    for(n in 1:n1){
    vector[L] p1_raw;
      for(l in 1:L){
        p1_raw[l] = beta1[l] * exp(normal_lpdf(y1[n] | mu[l], 1));
      }
      for(l in 1:L){
        p1[n, l] = p1_raw[l]/sum(p1_raw);
      }
  }
      for(n in 1:n2){
    vector[L] p2_raw;
      for(l in 1:L){
        p2_raw[l] = beta2[l] * exp(normal_lpdf(y2[n] | mu[l], 1));
      }
      for(l in 1:L){
        p2[n, l] = p2_raw[l]/sum(p2_raw);
      }
  }
}

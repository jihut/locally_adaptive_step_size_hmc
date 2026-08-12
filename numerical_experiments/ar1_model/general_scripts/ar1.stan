data {
  int d;
  real<lower=-1, upper=1> rho;
}

parameters {
  vector[d] x;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  // x[1] ~ normal(0, sqrt(1 / (1 - pow(rho, 2))));
  x[1] ~ normal(0, 1);
  // x[2:d] ~ normal(rho * x[1:(d - 1)], 1);
  x[2:d] ~ normal(rho * x[1:(d - 1)], sqrt(1 - pow(rho, 2)));
}

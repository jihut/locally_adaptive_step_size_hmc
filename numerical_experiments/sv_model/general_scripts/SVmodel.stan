
// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  vector[N] y;
  real mu0;
  real sd0;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real tRho;
  real logSigma;
  real x0;
  vector[N] z;
}


transformed parameters{

  real sigma = exp(logSigma);
  real rho = -1.0 + 2.0*exp(tRho)/(1.0+exp(tRho));
  
  vector[N] x; //(0...,T-1)
  vector[N] svsd;
  vector[N] mns;
  vector[N] sds;
  
  x[1] = x0;
  for(t in 2:N) x[t] = x[t-1] + sigma*z[t-1];
  svsd = exp(0.5*x); //(0,...,T-1)
  mns = rho*(z.*svsd);
  sds = sqrt(1.0-rho*rho)*svsd;
}



// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  // standard prior for sigma
  target += -10.0*logSigma - exp(-2.0*logSigma)/20.0;
  // uniform(-1,1) prior for rho
  target += -2.0*log_sum_exp(-0.5*tRho,0.5*tRho);
  // standard normal prior for standardized innovations
  z ~ normal(0,1);
  
  // initial period prior
  x0 ~ normal(mu0,sd0);
  
  // likelihood
  y ~ normal(mns,sds);
  
  
}


/*
  This Stan program conducts a non-linear regression of data that is presented
in histogram format. The non-linear models used are modified versions of the
Gibbs distribution. The normalization constants for these distributions have a
very complex closed form. As such, the functions have to be normalized
numerically.
*/

functions {
    real model_dist(real x,         // Function argument
                    real xc,        // Complement of function argument on the
                                    // domain (defined later)
                    real[] theta,   // Parameters
                    real[] x_r,     // Data (real)
                    int[] x_i) {    // Data (integer)
      real T = theta[1];            // Thermal portion temperature
      real r0 = theta[2];           // Thermal/epithermal cross-over
      real alpha = theta[3];        // Pareto Exponent
      real v;
      v = 1 / (1 + (x / r0)^2)^((1 + alpha) / 2) * exp(-r0/T * atan(x / r0));
      return v;
    }
    real model_int( real[] limits, real[] theta, data real[] x_r) {
      int x_i[0];
      real lowlim = limits[1];
      real uplim = limits[2];
      return integrate_1d(model_dist, lowlim, uplim, theta, x_r, x_i, 1.49e-8);
    }
}

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;             // number of histogram bins
  vector[N] x;                // histogram bin lower bound
  vector[N] y;                // histogram bin occupancy
}

transformed data {
  real x_r[0];
  vector[N] y_norm;
  y_norm = log(y / sum(y));   // log norm the bin occupancy into a probability
}

// The parameters accepted by the model.
parameters {
  real<lower = 0> sigma;      // prediction error scale
  real<lower = 0> T;          // temperature of thermal region
  real<lower = 0> r0;         // cross-over point
  real<lower = 0> alpha;      // pareto exponent of epithermal region
}

// Generate the predicted values from the parameter space
transformed parameters {
  vector[N] y_hat;
  for (n in 1:N)
    if (n < N)
      y_hat[n] = model_int({x[n], x[n + 1]}, {T, r0, alpha}, x_r);
  else
    y_hat[n] = model_int({x[n], positive_infinity()}, {T, r0, alpha}, x_r);
  y_hat = log(y_hat / sum(y_hat)); // normalization of estimated distribution
}

// The model to be estimated. We model the output
// 'y_norm' to be normally distributed with mean 'y_hat'
// and standard deviation 'sigma'.
model {
  sigma ~ gamma(2, 8e-6);
  T ~ gamma(4, 8e-5);
  r0 ~ gamma(4, 4e-5);
  alpha ~ gamma(2, 1);
  y_norm ~ normal(y_hat, sigma);
}

// Here we do the simulations from the posterior predictive distribution
generated quantities {
  real lpd;
  lpd = normal_lpdf(y_norm | y_hat, sigma);
  /*
    vector<lower = 0>[N] y_rep; // vector of same length as the data y
    for (n in 1:N)
      y_rep[n] = normal_rng(y_hat[n], sigma);
    */
}

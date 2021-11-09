functions {
  matrix designMatrix(
      int l, 
      vector log_mu,
      real rbf_variance,
      vector ml,
      int q
    ) {

    real h;
    matrix [q, l] X;
    h = (ml[2] - ml[1]) * rbf_variance;
    for (i in 1:q) {
      for (j in 1:l) {
        X[i, j] = 1;
      }
    }
    X[, 2] = log_mu;
    for (i in 1:(l - 2)) {
      vector[q] tmp;
      for (j in 1:q) {
        tmp[j] = pow(log_mu[j] - ml[i], 2);
      }
      X[, i + 2] = exp(-0.5 * tmp / pow(h, 2));
    }
    return(X);
  }
}

data {
  int q;
  int n;
  int p;
  int<lower=0> counts[q, n];
  vector[q] mu_mu;
  real smu;
  real astwo;
  real bstwo;
  int l;
  vector [n] size_factors;
  vector [l] mbeta;
  matrix[l, l] vbeta;
  real rbf_variance;
  vector[l - 2] ml;
  real eta;
}

transformed data {
  vector [n] log_size_factors = log(size_factors);
}

parameters {
  vector [q] log_mu;
  vector <lower=0> [q] delta;
  vector [l] beta;
  real <lower=0> stwo;
  vector <lower=0> [q] lambda;
}

transformed parameters {
  vector [q] fu = designMatrix(l, log_mu, rbf_variance, ml, q) * beta;
  vector <lower=0> [q] mu = exp(log_mu);
  vector [q] epsilon = log(delta) - fu;
}

model {
  stwo ~ inv_gamma(astwo, bstwo);
  beta ~ multi_normal(mbeta, stwo * vbeta);
  log_mu ~ normal(mu_mu, smu);
  lambda ~ gamma(eta / 2, eta / 2);
  delta ~ lognormal(fu, stwo ./ sqrt(lambda));
  for (j in 1:n) {
    counts[, j] ~ neg_binomial_2_log(
      log_size_factors[j] + log_mu,
      1 ./ delta
    );
  }
}

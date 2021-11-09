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
  int sq;
  int p;
  int<lower=0> counts[q, n];
  int<lower=0> spikes[sq, n];
  real as;
  real bs;
  real atheta;
  real btheta;
  matrix [n, p] batch_design;
  vector [n] aphi;
  vector[q] mu_mu;
  real smu;
  real astwo;
  real bstwo;
  int l;
  vector [l] mbeta;
  matrix[l, l] vbeta;
  real rbf_variance;
  vector[l - 2] ml;
  real eta;
  vector [sq] spike_levels;
}

parameters {
  vector [q] log_mu;
  vector <lower=0> [q] delta;
  simplex[n] tphi;
  real <lower=0> nu[n];
  vector <lower=0> [n] s;
  vector <lower=0> [p] theta;
  vector [l] beta;
  real <lower=0> stwo;
  vector <lower=0> [q] lambda;
}

transformed parameters {
  vector <lower=0> [n] theta_vector = batch_design * theta;
  vector <lower=0> [n] phi = tphi * n;
  vector <lower=0> [q] mu = exp(log_mu);
  vector [q] fu = designMatrix(l, log_mu, rbf_variance, ml, q) * beta;
  vector [q] epsilon = log(delta) - fu;
}


model {
  theta ~ gamma(atheta, btheta);
  tphi ~ dirichlet(aphi);
  stwo ~ inv_gamma(astwo, bstwo);
  beta ~ multi_normal(mbeta, stwo * vbeta);
  log_mu ~ normal(mu_mu, smu);
  lambda ~ gamma(eta / 2, eta / 2);
  delta ~ lognormal(fu, stwo ./ sqrt(lambda));
  s ~ gamma(as, bs);
  nu ~ gamma(1 ./ theta_vector, 1 ./ (s .* theta_vector));
  for (j in 1:n) {
    counts[, j] ~ neg_binomial_2_log(
      log(phi[j]) + log(nu[j]) + log_mu,
      1 ./ delta
    );
    spikes[, j] ~ poisson(nu[j] * spike_levels);
  }
}

data {
  int<lower=0> N;
  int<lower=0> p;
  matrix[N, p] X;
  vector[N] y;
}
parameters {
  real<lower=0> sigmaSq;
  vector[p] z; // standard normal
  vector<lower=0>[p] nu; // gamma(a, 1)
  vector<lower=0>[p] lambda; // gamma(b, 1)
  real<lower=0> a;
  real<lower=0> b;
  real<lower=0> phi;
}

transformed parameters {
  vector[p] beta;
  vector[p] local_scale;
  local_scale = sqrt(phi * nu ./ lambda);
  beta = z .* local_scale;
}

model {
  sigmaSq ~ inv_gamma(1.5, 0.5);
  a ~ gamma(1.5, 1);
  b ~ gamma(1.5, 1);
  phi ~ cauchy(0, 1);

  nu ~ gamma(a, 1);
  lambda ~ gamma(b, 1);
  z ~ std_normal();
  y ~ normal(X * beta, sqrt(sigmaSq));
}

functions {
  real[] de2(real t, real[] y, real[] p, real[] r, int[] i) {
    int nn = i[1];
    matrix[nn, 2] dY;
    real tau = p[1];
    vector[nn] a = to_vector(p[2:nn+1]);
    vector[nn] v = to_vector(y[:nn]);
    vector[nn] w = to_vector(y[nn+1:]);
    dY[:, 1] = (v - square(v).*v/3.0 + w) * tau;
    dY[:, 2] = (a - v) / tau;
    return to_array_1d(dY);
  }
}

data {
  int<lower=1> nn;
  int<lower=0> nc;
  int ci[nc];
  real ct[nc];
  real cmu[nc];
  real csd[nc];
}

transformed data {
  int ns = 2 * nn;
  real tlo = 0.5;
  real thi = 5;
  real ode_xr[0];
  int ode_xi[1] = {nn};
}

parameters {
  real<lower=tlo, upper=thi> tau;
  real a[nn];
  real<lower=-10, upper=10> y0[ns];
}

model { 
  real ys[nc, ns];
  real ode_par[1 + nn];
  a ~ normal(1.05, 0.1);
  tau ~ normal(3, 0.2) T[tlo, thi];
  ode_par[1] = tau;
  ode_par[2:nn+1] = a;
  ys = integrate_ode_rk45(de2, y0, 0.0, ct, ode_par, ode_xr, ode_xi); 
  for (i in 1:nc)
    ys[i, ci[i]] ~ normal(cmu[i], csd[i]);
}

generated quantities {
  real ys[100, ns];
  real ts[100];
  real ode_par[1 + nn];
  ode_par[1] = tau;
  ode_par[2:nn+1] = a;
  for (i in 1:100)
    ts[i] = i/10.0;
  ys = integrate_ode_rk45(de2, y0, 0.0, ts, ode_par, ode_xr, ode_xi); 
}

# vim: ft=stan sw=2 sts=2 et ai

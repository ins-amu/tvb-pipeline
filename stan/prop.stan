functions {
  vector diff_coupling(matrix fc, vector v) {
    int nn = rows(v);
    vector[nn] dc = rep_vector(0.0, nn);
    for (i in 1:nn)
      for (j in 1:nn)
        dc[i] = dc[i] + fc[i, j] * (v[j] - v[i]);
    return dc;
  }

  real[] de2(real t, real[] y, real[] p, real[] r, int[] i) {
    int nn = i[1];
    matrix[nn, 2] dY;
    real tau = p[1];
    real k = p[2];
    vector[nn] a = to_vector(p[3:nn+2]);
    matrix[nn, nn] fc = to_matrix(p[nn+2:nn+2+nn*nn], nn, nn);
    vector[nn] v = to_vector(y[:nn]);
    vector[nn] w = to_vector(y[nn+1:]);
    dY[:, 1] = (v - square(v).*v/3.0 + w) * tau;
    dY[:, 2] = (a - v + k * diff_coupling(fc, v)) / tau;
    print(dY);
    return to_array_1d(dY);
  }

  real[] ode_par(real tau, real k, real[] a, matrix fc) {
    int nn = size(a);
    real par[2 + nn + nn*nn];
    par[1] = tau;
    par[2] = k;
    par[3:nn+2] = a;
    par[nn+2:nn+2+nn*nn] = to_array_1d(fc);
    return par;
  }
}

data {
  int<lower=1> nn; // num nodes
  int<lower=0> nc; // num constraints
  int<lower=1, upper=nn*2> ci[nc]; // node svar idx / constraint
  real<lower=0.0> ct[nc]; // node time / constraint
  real cmu[nc]; // mu for svar ci at ct
  real<lower=0.0> csd[nc]; // sd for svar ci at ct
}

transformed data {
  int ns = 2 * nn;
  real tlo = 0.5;
  real thi = 5;
  real ode_xr[0];
  int ode_xi[1] = { nn };
}

parameters {
  real<lower=tlo, upper=thi> tau;
  real<lower=0.5, upper=1.5> a[nn];
  real<lower=0, upper=2> y0[ns];
  real<lower=0, upper=1> k;
  matrix<lower=0, upper=1>[nn, nn] fc;
}

transformed parameters {
  real par[2 + nn + nn*nn] = ode_par(tau, k, a, fc);
  real ys[nc, ns];
  ys = integrate_ode_rk45(de2, y0, 0.0, ct, par, ode_xr, ode_xi);
}

model { 
  a ~ normal(1.05, 0.1);
  tau ~ normal(3, 0.2) T[tlo, thi];
  //for (i in 1:nc)
  //  cmu[i] ~ normal(ys[i, ci[i]], csd[i]);
}

generated quantities {
  real ts[10];
  real ysf[nc, ns];
  ysf = integrate_ode_rk45(de2, y0, 0.0, ts, par, ode_xr, ode_xi);
}

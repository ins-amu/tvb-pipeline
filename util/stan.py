# coding: utf-8
# conda create -y -n pystan python=3 cython numpy nomkl scipy matplotlib jupyter pandas
# source activate pystan 
# pip install pystan seaborn

import numpy as np
import pystan
import time
import seaborn as sns
from util.pkl import write_pkl, read_pkl


def compile_model(stan_fname, pkl_fname):
    write_pkl(pystan.StanModel(file=stan_fname))


def optimize(model, data, n=1):
    return [model.optimizing(data=data) for _ in range(n)]


def optimize_files(model_pkl, data_pkl, out_pkl, n=1):
    n = int(n)
    model = read_pkl(model_pkl)
    data = read_pkl(data_pkl)
    write_pkl(out_pkl, optimize(model, data, n=n))


def sample(model, data, optim):
    fit = model.sampling(data=data, chains=len(optim), init=optim)
    est = fit.extract(permuted=True)
    return fit, est


def sample_files(model_pkl, data_pkl, optim_pkl, est_pkl):
    model = read_pkl(model_pkl)
    data = read_pkl(data_pkl)
    optim = read_pkl(optim_pkl)
    fit, est = sample(model, data, optim)
    write_pkl(est_pkl, est)


if __name__ == '__main__':
    import sys
    cmd = sys.argv[1]
    eval(cmd)(*sys.argv[2:])

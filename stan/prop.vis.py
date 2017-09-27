import sys
_, samp_pkl_fname, png_fname = sys.argv

import pylab as pl
from util.pkl import read_pkl

samp = read_pkl(samp_pkl_fname)
pl.figure()
pl.subplot(232)
pl.hist(samp['a'])

pl.subplot(231)
pl.hist(samp['tau'])

pl.show()
pl.savefig(png_fname)
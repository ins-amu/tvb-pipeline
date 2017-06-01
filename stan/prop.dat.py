import sys
import numpy as np
from util.pkl import write_pkl, read_pkl

data = {
    'ct': np.r_[4.7, 4.8],
    'cmu': np.r_[-1.5, -1.5],
    'csd': np.r_[0.3, 0.4],
    'ci': np.r_[1, 3],
    'nc': 2,
    'nn': 2,
}

pkl_fname = sys.argv[1]
write_pkl(pkl_fname, data)

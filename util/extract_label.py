import sys

import numpy as np
import nibabel as nib


infile = sys.argv[1]
label = int(sys.argv[2])
outfile = sys.argv[3]

mgz_in = nib.load(infile)

labels_in = mgz_in.get_data()
labels_out = np.zeros_like(labels_in)
labels_out[labels_in == label] = 1

mgz_out = nib.freesurfer.mghformat.MGHImage(labels_out, mgz_in.affine, mgz_in.header)
nib.save(mgz_out, outfile)

# recon make

This is a simplified pipeline using just make and a single Python module
for custom image processing algorithms.

Usage:

```bash
source activate

# make freesurfer recon, tractography and connectome
make SUBJECT=ac-p15 T1=data/t1/epi2009/ac-p15 DWI=data/dti/epi2009/ac-p15 fs-recon resamp-anat tck conn

# reconstruct electrodes
make SUBJECT=ac-p15 ELEC=data/ct/ac-p15.nii.gz elec_mode=ELEC labeled_elec
```

For rationale for using `make`, cf http://zmjones.com/make/

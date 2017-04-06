# README

This is a Makefile for building brain models for TVB.

Essential variables can be provided by environment variables or
arguments to `make`:

```bash
source activate
make T1=myT1.nii.gz DWI=dwi_dicom_folder fs-recon dwi-raw
```

*Required software*

- Python w/ NumPy, SciPy, NiBabel
- FreeSurfer 6
- FSL
- MRtrix3 (at least v0.3.15)
- OpenMEEG


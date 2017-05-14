# README

This is a Makefile for building brain models for TVB.

Essential variables can be provided by environment variables or
arguments to `make`:

```bash
source activate
make SUBJECT=cj elec_mode=CT
```
assuming you put data in `data/cj/T1.nii.gz`, `data/cj/DWI.mif`, etc. Otherwise specify
```
make ... CT=/path/to/CT.mgz
```

*Required*

- Python w/ NumPy, SciPy, NiBabel
- FreeSurfer 6
- FSL
- MRtrix3 (>= v0.3.15)

*providing an sEEG schema*

I was labeling the image, and then the schema would map
label values to names, but, better:

Provide a file with xyz locations in the original electrode
image (CT or T1) and electrode names.  Once the labeled
elec image is ready, xyz locs are mapped to nearest object
and corresponding name used to generate contact names.
Optionally with number to know how many contacts to generate.


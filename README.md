# tvb-make

This is a Makefile and supporting Python module for preprocessing structural data
for TVB brain models. The dataflow implemented can be seen in the following diagram:

![dag](dag.png)

## Usage

```bash
source activate

# make freesurfer recon, tractography and connectome
make SUBJECT=ac-p15 T1=data/t1/epi2009/ac-p15 DWI=data/dti/epi2009/ac-p15 fs-recon resamp-anat tck conn

# reconstruct electrodes
make SUBJECT=ac-p15 ELEC=data/ct/ac-p15.nii.gz elec_mode=ELEC labeled_elec
```

The `run-one.sh` script can be used with OAR to run many jobs, e.g.
```bash
oarsub -l nodes=1,walltime=24:00:00 --array-param-file params.txt ./run-one.sh
```

### Special Cases

#### JPEG encoded images

If you DICOM files are encoded with lossless JPEG compression, most
of the software used will fail to read them. You can have the pipeline
decompress those images prior to processing them by placing the
`.dcmdjpeg.dir` suffix on the DICOM directory. For example, if your
T1 DICOM files are in `data/t1`, you can specify
```
make T1=data/t1.dcmdjpeg.dir
```
and the files will be decompressed into the `data/t1.dcmdjpeg.dir`
directory prior to processing.

#### ADNI data

Diffusion data from the [ADNI](http://www.adni-info.org/) project
require stack the Nifti files and extract the gradient scheme from XML
files, which can be automated by renaming the DWI data directory with
the `.adni.dir` suffix and converting to Mrtrix image format via
```bash
mv dti_dir dti.adni.dir
make dti.mif
```
Alternatively, `dti.mif` can be provided as the `DWI` argument directly
and conversion is performed automatically,
```bash
mv dti_dir dti.adni.dir
make DWI=dti.mif T1=... fs-recon conn
```

## Dependencies

The present [Dockerfile](Dockerfile) should now be considered a complete,
operational specification of the pipeline's dependencies, but roughly
speaking, we require

- make
- Python w/ NumPy, SciPy, NiBabel, MNE
- FreeSurfer 6
- FSL
- MRtrix3 (>= v0.3.15)
- OpenMEEG
- DCMTK (optional, for decompressing JPEG DICOMs, see above)

It's likely easier to use our prebuilt Docker image.
[Install Docker](https://docs.docker.com/engine/installation/), pull the
image (`docker pull maedoc/tvb-make`) and run commands with the
Docker container (examples to come).

The main caveat is the Docker image doesn't do fancy graphics, so you'll
still want a native copy of Mrview, Freeview etc for visualization.

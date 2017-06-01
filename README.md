# tvb-make

This is a workflow for preparing and using TVB brain network models, comprised of
three main components

- process & dataflow in a [Makefile](Makefile)
- supporting Python utilities in a [util](util) module
- a [Dockerfile](docker/Dockerfile) & [image](https://hub.docker.com/r/maedoc/tvb-make/) with all dependencies

**Table of contents**

- [Functionality](#functionality)
- [Dependencies](#dependencies)
- [Usage](#usage)
  - [Targets](#targets)
  - [Docker](#docker)
  - [Marseille Cluster](#marseille-cluster)
- [Special Cases](#special-cases)
  - [JPEG encoded images](#jpeg-encoded-images)
  - [ADNI data](#adni-data)

## Functionality

- Automatic generation of surface & connectivity datasets usuable in TVB
- Forward models for MEG, EEG, sEEG with OpenMEEG
- Data fitting & parameter tuning via sensitivity analysis & Bayesian inversion

Many aspects are currently works in progress, but
the dataflow currently implemented can be seen in the following diagram:
![dag](dag.png)

## Dependencies

- make
- Python w/ NumPy, SciPy, NiBabel, MNE
- FreeSurfer 6
- FSL
- MRtrix3 (>= v0.3.15)
- OpenMEEG
- DCMTK (optional, for decompressing JPEG DICOMs, see above)

It's likely easier to use our prebuilt Docker image.
[Install Docker](https://docs.docker.com/engine/installation/), and
run commands with the Docker container (examples below).

The main caveat is the Docker image doesn't do fancy graphics, so you'll
still want a native copy of Mrview, Freeview etc for visualization.

## Usage

Basic usage requires invoked `make` with a subject name and your dataset,
```bash
make SUBJECT=tvb T1=data/t1 DWI=data/dwi fs-recon conn
```
where arguments are provided in `ARG=value` form, and outputs are given
as names like `fs-recon` to perform the FreeSurfer `recon-all -all`
reconstruction. See the following Targets section for a list of available
outputs.

### Targets

- `fs-recon` - FreeSurfer reconstruction, `recon-all -all ...`
- `resamp-anat` - Lower resolution cortical surfaces & annotations
- `conn` - Connectivity files

### Docker

Docker requires you to specify where your data is, so assuming you
have your data in `/home/me/data`, you might start with
```sh
docker run --rm -it -v /home/me/data:/data maedoc/tvb-make bash
```
and invoke commands as required.  Making this more user friendly
is still a work in progress.

### Marseille Cluster

Most dependencies are available on the Marseille cluster, and they
can be made available with
```bash
source activate
```

Additionally, the [run-one.sh](run-one.sh) script simplifies submitting
jobs to OAR, e.g.

The `run-one.sh` script can be used with OAR to run many jobs, e.g.
```bash
oarsub -l nodes=1,walltime=24:00:00 --array-param-file params.txt ./run-one.sh
```

## Special Cases

### JPEG encoded images

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

### ADNI data

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

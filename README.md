# pipeline

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

- `fs-recon`: FreeSurfer reconstruction. Consists mainly of running `recon-all -all`.
  Uses `T1`.

- `resamp-anat`: Lower resolution cortical surfaces & annotations
  Uses `T1`.

- `conn`: Connectivity matrices in text format.
  Uses `T1` and `DWI`.

- `tvb`: TVB zipfile, cortical and subcortical surfaces in TVB formats, region mappings.
  Uses `T1` and `DWI`.

- `elec`: Positions of the contacts of depth electrodes and gain matrices.
  Uses `T1`, `DWI`, `ELEC`, and either `ELEC_ENDPOINTS` or `ELEC_POS_GARDEL`.

- `seeg`: Conversion of SEEG recordings to FIF format and plotting the recordings.
  Uses `SEEGRECDIR`, `XLSX` and everything that `elec` uses.

- `ez`: Extraction of the epileptogenic zone from the patient Excel file.
  Uses `XLSX` and everything that `elec` uses.
 

_TODO_ more details & help on this

### Config file

As an added convenience, a file in Make format can be provided via the `CONFIG`
variable, with the desired values or even extra rules, overriding the defaults.
For example, the line
```bash
make SUBJECT=tvb T1=data/t1 fs-recon conn
```
could be replaced by a file `tvb.config.mk`
```make
SUBJECT := tvb
T1 := data/t1
.DEFAULT: fs-recon conn
```
and the invocation
```bash
make CONFIG=tvb.config.mk
```

### Docker

The `docker/run` script facilitates invoking the pipeline
in a virtual machine, so that no installation is required:
```bash
docker/run arguments...
```
The `data` folder in the current folder is available to the
container under the same name; place input data there and
provide corresponding paths to the pipeline. 
For example, if you use `/work/project1` as a working directory,
create `/work/project1/data`, place a T1 at `work/project1/data/T1.nii.gz`
and invoke as follows
```bash
~ $ cd /work/project1
/work/project1 $ /path/to/tvb-make/docker/run SUBJECT=tvb T1=data/T1.nii.gz fs-recon
...
```
The mapped directory can be customized with the `TVB_MAKE_DATA`
environment variable.

### Marseille Cluster

For quick introduction look at the basic [step-by-step tutorial](doc/TutorialCluster.md).

There are two options for running the pipeline on the cluster: non-interactive and interactive.
For running the full pipeline, non-interactive mode is recommended due to the large time requirements. 
For small updates and testing the interactive mode might be more suitable.

#### Non-interactive mode

In the non-interactive regime, you prepare the data and submit the job(s), and the scheduler takes cares of the execution.
The `cluster/run` script assists in running the pipeline on the Marseille cluster through two modes.
 First, invoke with typical arguments
```bash
<PIPELINE_DIR>/cluster/run SUBJECTS_DIR=fs SUBJECT=foo T1=data/T1.nii.gz fs-recon
```
for a single run in a single SLURM job. If you have many subjects, create a file `params.txt` with multiple lines of arguments, e.g.
```
SUBJECTS_DIR=fs SUBJECT=foo T1=data/T1.nii.gz fs-recon
SUBJECTS_DIR=fs SUBJECT=bar T1=data/T2.nii.gz fs-recon conn
SUBJECTS_DIR=fs SUBJECT=baz T1=data/T3.nii.gz conn
```
then
```
<PIPELINE_DIR>/cluster/run params.txt
```
Each line will result in the pipeline running a SLURM job for every line. You can comment out a line if you prepend it with a `#` sign,
```
 # SUBJECTS_DIR=fs SUBJECT=foo T1=data/T1.nii.gz fs-recon
```

NB You need to provide a custom, valid FreeSurfer `SUBJECTS_DIR`,
since the default directories on the cluster (`/soft/freesurfer*/subjects`)
are not writeable by users.

#### Interactive mode

First, request a computational node in the interactive mode
```
srun --pty bash
```
which should give you the interactive node if there is one available.

If you need to run the reconstruction and tractography in the interactive mode
(although that is discouraged), you need to request full node with enough memory:
```
srun -N 1 -n 1 --exclusive --mem=60G --pty bash
```

Then setup your working environment by loading the environment file,
```
source <PIPELINE_DIR>/cluster/env
```
and run `make` by hand:
```
make -f <PIPELINE_DIR>/Makefile  SUBJECTS_DIR=fs SUBJECT=foo T1=data/T1.nii.gz fs-recon
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

## Stan support

Generic support for Stan models is implemented in 
[`make/Stan.mk`](make/Stan.mk) with the following conventions:
for each Stan model, there are three files to provide in this repository:

- `stan/{model_name}.stan` - the Stan code for the model
- `stan/{model_name}.dat.py` - Python script to generate input data
- `stan/{model_name}.vis.py` - Python script to visualize results

which generate or use files in the `stan` subfolder of the subjects' folder
(`$(sd)` in the following):

- `$(sd)/stan/{model_name}.stan.pkl` - compiled Stan model in PyStan pickle format
- `$(sd)/stan/{model_name}.dat.pkl` - input data generated by `stan/{model_name}.dat.py`
- `$(sd)/stan/{model_name}.opt.pkl` - posterior mode found in initial optimization
- `$(sd)/stan/{model_name}.samp.pkl` - posterior samples found during fit
- `$(sd)/stan/{model_name}.png` - visualization produced by `stan/{model_name}.vis.py`

See the [`stan`](stan) folder for an example, to be completed.

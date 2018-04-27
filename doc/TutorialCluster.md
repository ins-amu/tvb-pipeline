


# Running the pipeline on INS cluster


Brief tutorial on how to run the reconstruction pipeline on the INS cluster.

## Setting up the environment (do once)

Clone the pipeline repository. Note that you need your SSH key generated on the cluster to be added
in the Gitlab interface to be able to clone the repo from the cluster.
```
cd ~/soft       # Or wherever you want the code to be
git clone git@gitlab.thevirtualbrain.org:tvb/pipeline.git
```

## Preparing the data

By default, the pipeline uses following data structure:
```
data/SUBJECT-1/
     SUBJECT-2/
     SUBJECT-3/
     ...
fs/SUBJECT-1/
   SUBJECT-2/
   SUBJECT-3/
   ...
```
where `data/` contains the raw data, and `fs/` contains the processed results. You should prepare
the contents of the `data/` directory; the contents of the `fs/` directory is filled by the pipeline.
If needed, you can change the the names of the raw data directory by setting the make variable `DATA`
(see below on how to), and the `fs/` directory by the variable `SUBJECTS_DIR`.


First, let's create this main directory structure in `~/reconstruction`
```
cd
mkdir -p reconstruction/data reconstruction/fs
```


Then start with a single patient `SUBJECT-1`. 
For the basic TVB dataset, you need at least T1 and DWI scans. Then place the T1 and DWI scans under `t1/`
and `dwi/` directories under the patient directory:
```
data/SUBJECT-1/dwi/
               t1/
```


## Running the pipeline

You need to create a file containing the subject specification. In your working directory, create
a file `params.txt` and insert a single line inside:
```
SUBJECT=SUBJECT-1 T1=data/SUBJECT-1/t1/t1.nii.gz DWI=data/SUBJECT-1/dwi/dwi.nii.gz tvb
```
Now what this means? Every line specifies a single subject, so here we have specified a single
subject named `SUBJECT-1`. 
By setting the variables `T1` and `DWI` we have told the pipeline where the raw data are.
If you have the raw data in DICOM format (many .dcm files in `t1` and `dwi` 
directories) and not single files, simply point to the directories: 
```
SUBJECT=SUBJECT-1 T1=data/SUBJECT-1/t1/ DWI=data/SUBJECT-1/dwi/ tvb
```

Last keyword on this line is the *target* of the pipeline. In this case, `tvb` stands for the 
TVB data set with connectomes and surfaces. Other targets that may be useful are `fs-recon` for
the FreeSurfer reconstruction only, `elec` for the depth electrode positions, or `seeg` for SEEG
recordings.


The pipeline job is submitted simply by running the following command:

```
~/soft/pipeline/cluster/run params.txt
```
The output of the command should show that for every line in the `params.txt` (not commented out) 
a Slurm job was submitted.

## Examining the results

The status of the Slurm jobs on the cluster can be checked by
```
squeue -u USERNAME
```
If you have submitted a job and it is finished after only a few seconds, something probably
went wrong. Have a look at the logs in `fs/_logs/SUBJECT-1.*.stdout` and
 `fs/_logs/SUBJECT-1.*.stderr`.

After the job has ended, examine the logs and the created subject directory
`fs/SUBJECT-1/`, especially the `tvb` subdirectory, where your desired data (TVB zipfiles,
 surfaces, region mappings) should be.







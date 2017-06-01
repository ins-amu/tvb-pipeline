#!/bin/bash

# Preprocessing of ADNI imaging files
# Usage:
# > export DTI_DIR=...
# > export DTI_MERGED_FILE=...
# > ./preproc_adni.sh

fslmerge -t $DTI_DIR/$DTI_MERGED_FILE $DTI_DIR/*.nii
python -m util.xml2bvalsbvecs $DTI_DIR/bvals.txt $DTI_DIR/bvecs.txt $DTI_DIR/*.xml

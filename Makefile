# Issue tracker ;) {{{ 
# TODO bvecs bvals
# TODO redo repo and make public
# TODO use include $(var) to split up
# TODO move activate here
# TODO redirect dwi & seeg to separate log files
# TODO use an isrunning file to avoid running dwi twice simultaneously
# TODO 2d parsweep simple stability subcort & cortical coupling scaling to recalibrate matrix
# TODO target OAR directly, stdout/stderr to subject directories
# TODO work with / inside of docker image
# }}}

sval ?= pial
pe_dir ?= AP
regopts ?= -cost mutualinfo -dof 12 -searchrz -180 180 -searchry -180 180  -searchrx -180 180 
ntrack ?= 15M
ct_thresh ?= 1000
nthread ?= 8
parc ?= aparc.a2009s
aa ?= aparc+aseg
resamp_target ?= fsaverage5
elec_mode ?= CT

# default data layout {{{
data ?= data
T1 ?= $(data)/$(SUBJECT)/T1.nii.gz
DWI ?= $(data)/$(SUBJECT)/DWI.mif
CT ?= $(data)/$(SUBJECT)/CT.nii.gz
ELEC ?= $(data)/$(SUBJECT)/ELEC.nii.gz
# }}}

hemi = lh rh
sd = $(SUBJECTS_DIR)/$(SUBJECT)
rtd = $(SUBJECTS_DIR)/$(resamp_target)
fs_done = $(sd)/mri/$(aa).mgz
here := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))

.PHONY: fs-recon resamp-anat dwi seeg clean mrinfo # {{{
default:
	echo "please read Makefile to use correctly"
fs-recon: $(fs_done) $(sd)/mri/$(aa).xyz
resamp-anat: $(sd)/surf/lh.$(sval).$(resamp_target) $(sd)/surf/rh.$(sval).$(resamp_target)
tck: $(sd)/dwi/100k.tck
conn: $(sd)/dwi/counts.txt $(sd)/dwi/lengths.txt
labeled_elec: $(sd)/seeg/labeled_$(elec_mode).nii.gz
seeg: $(sd)/seeg/seeg.xyz $(sd)/seeg/gain.mat

mrinfo:
	echo $(SUBJECT)
	mrinfo $(T1)
	mrinfo $(DWI)

clean:
	rm -rf $(sd) $(SUBJECT)-T1_raw.nii.gz

nothing:
	echo "doing nothing for subject $(SUBJECT) per request"

# }}}

# FreeSurfer reconstruction & downsampling# {{{

$(sd)/mri/orig/001.mgz: $(T1)
	mrconvert $(T1) $(SUBJECT)-T1_raw.nii.gz
	recon-all -s $(SUBJECT) -i $(SUBJECT)-T1_raw.nii.gz
	rm $(SUBJECT)-T1_raw.nii.gz

$(fs_done): $(sd)/mri/orig/001.mgz
	recon-all -s $(SUBJECT) -all -parallel -openmp $(nthread)

$(rtd):
	cp -r $(FREESURFER_HOME)/subjects/$(resamp_target) $(SUBJECTS_DIR)/

$(sd)/surf/%.$(sval).$(resamp_target): $(rtd) $(fs_done)
	mri_surf2surf \
		--srcsubject $(SUBJECT) \
		--trgsubject $(resamp_target) \
		--hemi $* \
		--sval-xyz $(sval) \
		--tval $(sval).$(SUBJECT) \
		--tval-xyz $(sd)/mri/T1.mgz
	cp $(rtd)/surf/$*.$(sval).$(SUBJECT) \
	    $(sd)/surf/$*.$(sval).$(resamp_target)
	mri_surf2surf \
		--srcsubject $(SUBJECT) \
		--trgsubject $(resamp_target) \
		--hemi $* \
		--sval-annot $(sd)/label/$*.$(parc).annot \
		--tval $(sd)/label/$*.$(parc).annot.$(resamp_target)

$(sd)/mri/$(aa).xyz: $(fs_done)
	./reconutil label_volume_centers $(sd)/mri/$(aa).mgz $@

# }}}

# Diffusion processing# {{{
$(sd)/dwi/t2d.mat: $(fs_done) $(sd)/dwi/bzero.nii.gz $(sd)/mri/T1.RAS.nii.gz
	flirt -ref $(sd)/dwi/bzero.nii.gz \
	    -in $(sd)/mri/T1.RAS.nii.gz \
	    -omat $(sd)/dwi/t2d.mat \
	    -out $(sd)/dwi/T1_in_bzero.nii.gz \
	    $(regopts)

$(sd)/dwi/aparc_aseg.nii.gz: $(sd)/dwi/t2d.mat $(sd)/mri/$(aa).RAS.RO.nii.gz
	flirt -applyxfm -interp nearestneighbour \
	    -in $(sd)/mri/$(aa).RAS.RO.nii.gz \
	    -ref $(sd)/dwi/bzero.nii.gz \
	    -init $(sd)/dwi/t2d.mat -out $@

# 001.mgz: wait for subject folder to exist before proceeding
$(sd)/dwi/raw.mif: $(DWI) $(sd)/mri/orig/001.mgz
	mkdir -p $(sd)/dwi
	mrconvert -force $< $@

$(sd)/dwi/preproc.mif: $(sd)/dwi/raw.mif
	# dwipreproc -force -rpe_none $(pe_dir) $< $@
	cp $< $@

$(sd)/dwi/bzero.mif: $(sd)/dwi/preproc.mif
	dwiextract -force -bzero $< $@

$(sd)/dwi/response.txt: $(sd)/dwi/preproc.mif
	dwi2response -nthreads $(nthread) -force tournier $< $@

$(sd)/dwi/mask.mif: $(sd)/dwi/preproc.mif
	dwi2mask -force $< $@

$(sd)/dwi/fod.mif: $(sd)/dwi/preproc.mif $(sd)/dwi/response.txt
	dwi2fod csd $^ $@ -force -nthreads $(nthread)

$(sd)/dwi/100k.tck: $(sd)/dwi/all.tck
	tckedit -number 100K $< $@

$(sd)/dwi/all.tck: $(sd)/dwi/fod.mif $(sd)/dwi/mask.mif
	tckgen $< $@ \
	    -mask $(sd)/dwi/mask.mif \
	    -seed_image $(sd)/dwi/mask.mif \
	    -number $(ntrack) -force -nthreads $(nthread)

$(sd)/dwi/label.mif: $(sd)/dwi/aparc_aseg.nii.gz
	labelconvert $< \
	    $(FREESURFER_HOME)/FreeSurferColorLUT.txt \
	    $(MRT3)/src/connectome/tables/fs_default.txt \
	    $@ -force

$(sd)/dwi/counts.txt: $(sd)/dwi/all.tck $(sd)/dwi/label.mif 
	tck2connectome $^ $@ -force -nthreads $(nthread)

$(sd)/dwi/lengths.txt: $(sd)/dwi/all.tck $(sd)/dwi/label.mif 
	tck2connectome $^ $@ -scale_length -stat_edge mean -force -nthreads $(nthread)

# }}}

# CT sEEG {{{

$(sd)/seeg: $(sd)/mri/orig/001.mgz
	mkdir -p $(sd)/seeg

$(sd)/seeg/CT_in_T1.nii.gz: $(CT) $(fs_done) $(sd)/mri/T1.RAS.nii.gz $(sd)/seeg
	mri_convert $< $(sd)/seeg/CT.nii.gz --out_orientation RAS
	flirt -in $(sd)/seeg/CT.nii.gz -ref $(sd)/mri/T1.RAS.nii.gz \
	    -omat $(sd)/seeg/CT_to_T1.mat \
	    -out $(sd)/seeg/CT_in_T1.nii.gz \
	    $(regopts)

$(sd)/seeg/mask.nii.gz: $(sd)/mri/brain.RAS.nii.gz $(sd)/seeg
	mri_binarize --i $< --o $@ --min 10 --erode 4

$(sd)/seeg/masked_CT.nii.gz: $(sd)/seeg/CT_in_T1.nii.gz $(sd)/seeg/mask.nii.gz
	mri_binarize --i $< --o $@ --min $(ct_thresh) \
	    --mask $(sd)/seeg/mask.nii.gz

$(sd)/seeg/dilated_CT.nii.gz: $(sd)/seeg/masked_CT.nii.gz
	mri_binarize --i $< --o $@ --min 0.5 --dilate 2 --erode 1

$(sd)/seeg/labeled_CT.nii.gz: $(sd)/seeg/masked_CT.nii.gz $(sd)/seeg/dilated_CT.nii.gz
	./reconutil label_with_dilation $^ $@

# case of T1 w/ elec, no dilation
$(sd)/seeg/ELEC_in_T1.nii.gz: $(ELEC) $(fs_done) $(sd)/mri/T1.RAS.nii.gz $(sd)/seeg
	mri_convert $< $(sd)/seeg/ELEC.nii.gz --out_orientation RAS
	flirt -in $(sd)/seeg/ELEC.nii.gz -ref $(sd)/mri/T1.RAS.nii.gz \
	    -omat $(sd)/seeg/ELEC_to_T1.mat \
	    -out $(sd)/seeg/ELEC_in_T1.nii.gz \
	    $(regopts)

$(sd)/seeg/masked_ELEC.nii.gz: $(sd)/seeg/ELEC_in_T1.nii.gz $(sd)/seeg/mask.nii.gz
	mri_binarize --i $< --o $@ --min $(elec_min) --max $(elec_max) --mask $(sd)/seeg/mask.nii.gz

$(sd)/seeg/labeled_ELEC.nii.gz: $(sd)/seeg/masked_ELEC.nii.gz
	./reconutil label_objects $< $@ $(nthread)

# user must create schema.txt
$(sd)/seeg/seeg.xyz: $(sd)/seeg/labeled_$(elec_mode).nii.gz $(sd)/seeg/schema.txt
	./reconutil gen_seeg_xyz $^ $@

$(sd)/seeg/gain.mat: $(sd)/seeg/seeg.xyz $(sd)/mri/$(aa).xyz
	./reconutil seeg_gain $^ $@

# }}}

# Conversion rules# {{{
%.RAS.nii.gz: %.mgz
	mri_convert -rt nearest --out_orientation RAS $< $@

%.RO.nii.gz: %.nii.gz
	fslreorient2std $< $@

%.nii.gz: %.mif
	mrconvert -force $< $@

%.mif: %.nii.gz
	mrconvert -force $< $@
# }}}

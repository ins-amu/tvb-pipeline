# default parameters {{{
sval ?= pial
pe_dir ?= AP
regopts ?= -cost mutualinfo -dof 12 -searchrz -180 180 -searchry -180 180  -searchrx -180 180
ntrack ?= 15M
ct_thresh ?= 1000
nthread ?= 8
parc ?= aparc.a2009s
aa ?= aparc+aseg
resamp_target ?= fsaverage5
elec_mode ?= ELEC
lut_fs := $(FREESURFER_HOME)/FreeSurferColorLUT.txt
lut_target ?= $(shell find $(MRT3) -name fs_default.txt | head -n 1)
ifneq ($(and $(BVECS),$(BVALS)),)
    raw_mif_convert_flags := -fslgrad $(BVECS) $(BVALS)
endif
# }}}

# default data layout {{{
data ?= data
T1 ?= $(data)/$(SUBJECT)/t1/
DWI ?= $(data)/$(SUBJECT)/dwi/
CT ?= $(data)/$(SUBJECT)/elec/CT.nii.gz
ELEC ?= $(data)/$(SUBJECT)/elec/elec_ct.nii.gz
ELEC_POS_GARDEL ?= $(data)/$(SUBJECT)/elec/pos_vox.txt
SEEGRECDIR ?= $(data)/$(SUBJECT)/seeg
XLSX ?= $(data)/patient.xlsx
# }}}

# misc util {{{
export SUBJECTS_DIR
hemi = lh rh
sd = $(SUBJECTS_DIR)/$(SUBJECT)
rtd = $(SUBJECTS_DIR)/$(resamp_target)
fs_done = $(sd)/mri/$(aa).mgz
here := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
export PYTHONPATH:=$(here):$(PYTHONPATH)
# }}}

.PHONY: fs-recon resamp-anat dwi seeg tvb clean mrinfo # {{{
default:
	echo "please read Makefile to use correctly"
fs-recon: $(fs_done) $(sd)/mri/$(aa).xyz
resamp-anat: $(sd)/surf/lh.$(sval).$(resamp_target) $(sd)/surf/rh.$(sval).$(resamp_target)
tck: $(sd)/dwi/100k.tck
conn: $(sd)/dwi/counts.txt $(sd)/dwi/lengths.txt
labeled_elec: $(sd)/elec/labeled_$(elec_mode).nii.gz
elec: $(sd)/elec/seeg.xyz $(sd)/elec/gain_dipole_no-subcort.mat $(sd)/elec/gain_inv-square.mat $(sd)/elec/img
tvb: $(sd)/tvb/connectivity.zip $(sd)/tvb/img/connectivity.png

mrinfo:
	echo $(SUBJECT)
	mrinfo $(T1)
	mrinfo $(DWI)

clean:
	rm -rf $(sd) $(SUBJECT)-T1_raw.nii.gz

nothing:
	echo "doing nothing for subject $(SUBJECT) per request"

dag.png: Makefile
	mkdir dag
	touch dag/{t1,dwi}
	make -Bnd SUBJECTS_DIR=subjects SUBJECT=dag T1=t1 DWI=dwi fs-recon conn | make2graph | dot -Tpng -o dag.png
	rm -r dag

docker:
	cd docker && sudo docker build -t maedoc/tvb-make .

tmux:
	tmux

notebook:
	jupyter notebook --allow-root --ip=0.0.0.0

# }}}

# rules {{{
include $(here)/make/FreeSurfer.mk
include $(here)/make/Diffusion.mk
include $(here)/make/StereoEEG.mk
include $(here)/make/Stan.mk
include $(here)/make/Utilities.mk
#}}}

# user provided config file
-include $(CONFIG)

# vim: foldmethod=marker

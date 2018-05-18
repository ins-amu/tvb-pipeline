# default parameters {{{
pe_dir ?= AP
regopts ?= -cost mutualinfo -dof 12 -searchrz -180 180 -searchry -180 180  -searchrx -180 180
ntrack ?= 15M
ct_thresh ?= 1000
nthread ?= 8
resamp_target ?= fsaverage5
resamp_parc ?= aparc.a2009s
resamp_sval ?= pial
elec_mode ?= CT
lut_fs := $(FREESURFER_HOME)/FreeSurferColorLUT.txt
ifneq ($(and $(BVECS),$(BVALS)),)
    raw_mif_convert_flags := -fslgrad $(BVECS) $(BVALS)
endif
# }}}

# default data layout {{{
DATA ?= data

# Shorthand for the Data Directory for use in the params file
DD = $(DATA)/$(SUBJECT)

T1 ?= $(DATA)/$(SUBJECT)/t1/
DWI ?= $(DATA)/$(SUBJECT)/dwi/
ELEC ?= $(DATA)/$(SUBJECT)/elec/elec.nii.gz
SEEGRECDIR ?= $(DATA)/$(SUBJECT)/seeg
XLSX ?= $(DATA)/$(SUBJECT)/patient.xlsx

# Only one of the following should be present
ELEC_POS_GARDEL ?= $(DATA)/$(SUBJECT)/elec/pos_vox.txt
ELEC_ENDPOINTS ?= $(DATA)/$(SUBJECT)/elec/elec_endpoints.txt
ELEC_LABEL_SCHEMA ?= $(DATA)/$(SUBJECT)/elec/schema.txt
# }}}

# misc util {{{
export SUBJECTS_DIR
sd = $(SUBJECTS_DIR)/$(SUBJECT)
rtd = $(SUBJECTS_DIR)/$(resamp_target)
fs_done = $(sd)/mri/aparc+aseg.mgz
here := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
export PYTHONPATH:=$(here):$(PYTHONPATH)
# }}}

.PHONY: fs_recon resamp-anat dwi seeg tvb clean mrinfo # {{{
default:
	echo "please read Makefile to use correctly"
fs-recon: $(fs_done)
resamp-anat: $(sd)/surf/lh.$(resamp_sval).$(resamp_target) $(sd)/surf/rh.$(resamp_sval).$(resamp_target)
tck: $(sd)/dwi/100k.tck
conn: $(sd)/dwi/counts.dk.txt        $(sd)/dwi/lengths.dk.txt  \
	  $(sd)/dwi/counts.destrieux.txt $(sd)/dwi/lengths.destrieux.txt
conn_label: $(sd)/dwi/label_in_T1.dk.nii.gz $(sd)/dwi/label_in_T1.destrieux.nii.gz
labeled_elec: $(sd)/elec/labeled_elec-$(elec_mode).nii.gz
elec: $(sd)/elec/seeg.xyz $(sd)/elec/img \
      $(sd)/elec/elec.dk.png $(sd)/elec/elec.destrieux.png \
      $(sd)/elec/gain_dipole_no-subcort.dk.txt $(sd)/elec/gain_dipole_no-subcort.destrieux.txt \
	  $(sd)/elec/gain_inv-square.dk.txt        $(sd)/elec/gain_inv-square.destrieux.txt
tvb: $(sd)/tvb/connectivity.dk.zip         $(sd)/tvb/img/connectivity.dk.png \
     $(sd)/tvb/connectivity.destrieux.zip  $(sd)/tvb/img/connectivity.destrieux.png
ez: $(sd)/tvb/ez_hypothesis.dk.txt         $(sd)/tvb/img/ez_hypothesis.dk.png \
	$(sd)/tvb/ez_hypothesis.destrieux.txt  $(sd)/tvb/img/ez_hypothesis.destrieux.png

mrinfo:
	echo $(SUBJECT)
	mrinfo $(T1)
	mrinfo $(DWI)

clean:
	rm -rf $(sd) $(SUBJECT)-T1_raw.nii.gz

nothing:
	echo "doing nothing for subject $(SUBJECT) per request"

dag.png: Makefile make/FreeSurfer.mk make/Diffusion.mk make/StereoEEG.mk make/Stan.mk make/Utilities.mk
	mkdir -p dag/t1 dag/dwi dag/elec dag/seeg
	touch dag/elec/pos_vox.txt dag/elec/elec.nii.gz dag/patient.xlsx
	make -Bnd SUBJECTS_DIR=. SUBJECT=dag DATA=. \
		dag/mri/aparc+aseg.dk.mgz dag/mri/aparc+aseg.destrieux.mgz \
		dag/dwi/t2d.mat \
		dag/mri/aparc+aseg.dk.RAS.nii.gz dag/mri/aparc+aseg.destrieux.RAS.nii.gz \
		dag/mri/aparc+aseg.dk.RAS.RO.nii.gz dag/mri/aparc+aseg.destrieux.RAS.RO.nii.gz \
		dag/dwi/aparc+aseg.dk.nii.gz dag/dwi/aparc+aseg.destrieux.nii.gz \
		dag/dwi/label.dk.mif dag/dwi/label.destrieux.mif \
		dag/dwi/all.tck \
		dag/dwi/triu_counts.dk.txt dag/dwi/triu_counts.destrieux.txt \
		resamp-anat conn tvb elec seeg ez \
		| make2graph  \
		| egrep -v "Makefile|make/.*\.mk" \
		| sed '1 agraph [pad="0.5", nodesep="0.2", ranksep="1.7"];' \
		| dot -Tpng -o dag.png
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

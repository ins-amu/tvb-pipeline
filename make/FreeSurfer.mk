# Rules for FreeSurfer reconstruction & downsampling

# ensure user provided SUBJECTS_DIR propagates to commands
# as an environment variable
export SUBJECTS_DIR

# most other parts of the pipeline are organized into topic
# folders inside the subject folder, so we need FreeSurfer to
# setup that folder first; it just requires importing the
# T1. Other rules depend on the 001.mgz, so that they can
# run in parallel to the recon-all -all below
$(sd)/mri/orig/001.mgz: $(T1)
	tmp=$$(mktemp -d) && mrconvert $(T1) $$tmp/T1.nii.gz && recon-all -s $(SUBJECT) -i $$tmp/T1.nii.gz  && rm -r $$tmp
# run the full reconstruction w/ parallelism
$(fs_done): $(sd)/mri/orig/001.mgz
	recon-all -s $(SUBJECT) -all -parallel -openmp $(nthread)

$(sd)/mri/T1.RAS.nii.gz: $(fs_done)
	mri_convert -rt nearest --out_orientation RAS $(sd)/mri/T1.mgz $@


# this ensures we have a writeable copy of the resample target data
$(rtd):
	cp -r $(FREESURFER_HOME)/subjects/$(resamp_target) $(SUBJECTS_DIR)/

# resample anatomy, usually at lower resolution
$(sd)/surf/%.$(resamp_sval).$(resamp_target): $(rtd) $(fs_done)
	mri_surf2surf \
		--srcsubject $(SUBJECT) \
		--trgsubject $(resamp_target) \
		--hemi $* \
		--sval-xyz $(resamp_sval) \
		--tval $(resamp_sval).$(SUBJECT) \
		--tval-xyz $(sd)/mri/T1.mgz
	cp $(rtd)/surf/$*.$(resamp_sval).$(SUBJECT) \
	    $(sd)/surf/$*.$(resamp_sval).$(resamp_target)
	mri_surf2surf \
		--srcsubject $(SUBJECT) \
		--trgsubject $(resamp_target) \
		--hemi $* \
		--sval-annot $(sd)/label/$*.$(resamp_parc).annot \
		--tval $(sd)/label/$*.$(resamp_parc).annot.$(resamp_target)

# generate centers
$(sd)/mri/$(aa).xyz: $(fs_done)
	python -m util.util label_volume_centers $(sd)/mri/$(aa).mgz $@

$(sd)/mri/aparc+aseg.dk.mgz: $(fs_done)
	ln -sf ./aparc+aseg.mgz $@

$(sd)/mri/aparc+aseg.destrieux.mgz: $(fs_done)
	ln -sf ./aparc.a2009s+aseg.mgz $@

$(sd)/label/%.aparc.dk.annot: $(fs_done)
	cp $(sd)/label/$*.aparc.annot $@

$(sd)/label/%.aparc.destrieux.annot: $(fs_done)
	cp $(sd)/label/$*.aparc.a2009s.annot $@

# VEP parcellation/segmentation

$(sd)/mri/aparc+aseg.vep.mgz: $(sd)/mri/aparc+aseg.destrieux.mgz
	python -m util.convert_to_vep_parc convert_to_vep_parc   \
        $(sd)/mri/aparc+aseg.destrieux.mgz                   \
	    $(lut_fs) $(here)/util/data/VepAtlasRules.txt  $@

$(sd)/mri/aparc.vep.mgz: $(sd)/mri/aparc+aseg.vep.mgz
	python -m util.convert_to_vep_parc aparcaseg_to_aparc                     \
	    $(sd)/mri/aparc+aseg.vep.mgz $(here)/util/data/VepAparcColorLut.txt   \
	    $(here)/util/data/VepRegions.txt                                      \
	    $(sd)/mri/aparc.vep.mgz $(sd)/mri/aparc.lut.txt


$(sd)/label/%.aparc.vep.annot: $(sd)/mri/aparc.vep.mgz
	mris_sample_parc -ct $(here)/util/data/VepAparcColorLut.txt \
       $(SUBJECT) $* aparc.vep.mgz $*.aparc.vep.annot


# generate subcortical region volume bounding surfaces
$(sd)/aseg2srf/%: $(sd)/mri/aparc+aseg.%.mgz
	export UTILDIR=$(here)/util; \
	$(here)/util/aseg2srf -s $(SUBJECT) -l $(here)/util/data/subcort.$*.txt \
	    -a aparc+aseg.$*.mgz -o aseg2srf/$*

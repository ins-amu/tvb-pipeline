
# Conversion rules {{{
%.RAS.nii.gz: %.mgz
	mri_convert -rt nearest --out_orientation RAS $< $@

%.RO.nii.gz: %.nii.gz
	fslreorient2std $< $@

%.nii.gz: %.mif
	mrconvert -force $< $@

%.mif: %.nii.gz
	mrconvert -force $< $@

# e.g. make T1=/foo/bar.dcmdjpeg.dir ...
%.dcmdjpeg.dir: %
	mkdir -p $@
	cd $<; for img in *; do dcmdjpeg $$img ../$*.dcmdjpeg.dir/$$img; done

# put nii files in foo.adni.dir/; make DWI=foo.mif ...
%.mif: %.adni.dir
	fslmerge -t $*.nii.gz $</*.nii
	python -m util.xml2bvalsbvecs $</bvals.txt $</bvecs.txt $</*.xml
	mrconvert -fslgrad $</bvecs.txt $</bvals.txt $*.nii.gz $@

# }}}

# TVB compatible files {{{
$(sd)/tvb/connectivity.zip: $(fs_done) $(sd)/dwi/triu_counts.txt $(sd)/dwi/triu_lengths.txt $(sd)/aseg2srf
	mkdir -p $(sd)/tvb
	mris_convert $(sd)/surf/lh.pial $(sd)/surf/lh.pial.asc
	mris_convert $(sd)/surf/rh.pial $(sd)/surf/rh.pial.asc
	python -m util.create_tvb_dataset $(sd) \
	    $(lut_fs) $(lut_target) \
		$(sd)/dwi/triu_counts.txt $(sd)/dwi/triu_lengths.txt \
	    $(sd)/tvb/connectivity.zip $(sd)/tvb

$(sd)/tvb/img/connectivity.png: $(sd)/tvb/connectivity.zip
	mkdir -p $(sd)/tvb/img
	python -m util.plot plot_connectivity $< $@

# }}}

# Diffusion processing

# set this to -number for old mrtrix3
number_flag = -select

dwi_log := >> $(sd)/dwi/log.txt 2>&1

# 001.mgz: wait for subject folder to exist before proceeding
$(sd)/dwi/raw.mif: $(DWI) $(sd)/mri/orig/001.mgz
	mkdir -p $(sd)/dwi
	mrconvert $(raw_mif_convert_flags) -force $< $@ $(dwi_log)

# register T1 to DWI image
$(sd)/dwi/t2d.mat: $(fs_done) $(sd)/dwi/bzero.nii.gz $(sd)/mri/T1.RAS.RO.nii.gz
	flirt -ref $(sd)/dwi/bzero.nii.gz \
	    -in $(sd)/mri/T1.RAS.RO.nii.gz \
	    -omat $(sd)/dwi/t2d.mat \
	    -out $(sd)/dwi/T1_in_bzero.nii.gz \
	    $(regopts) $(dwi_log)

# move label volume to DWI space
$(sd)/dwi/aparc+aseg.%.nii.gz: $(sd)/mri/aparc+aseg.%.RAS.RO.nii.gz $(sd)/dwi/t2d.mat
	flirt -applyxfm -interp nearestneighbour \
	    -in $(sd)/mri/aparc+aseg.$*.RAS.RO.nii.gz  \
	    -ref $(sd)/dwi/bzero.nii.gz \
	    -init $(sd)/dwi/t2d.mat -out $@ $(dwi_log)

# preprocess DWI (TODO)
$(sd)/dwi/preproc.mif: $(sd)/dwi/raw.mif
	# dwipreproc -force -rpe_none $(pe_dir) $< $@
	# TODO reconutil func to invoke dwipreproc correctly
	cp $< $@ $(dwi_log)

# extract b0 volume
$(sd)/dwi/bzero.mif: $(sd)/dwi/preproc.mif
	dwiextract -force -bzero $< $@ $(dwi_log)

# estimate DWI response function
$(sd)/dwi/response.txt: $(sd)/dwi/preproc.mif
	dwi2response -nthreads $(nthread) -force tournier $< $@ $(dwi_log)

# extract brain mask volume
$(sd)/dwi/mask.mif: $(sd)/dwi/preproc.mif
	dwi2mask -force $< $@ $(dwi_log)

# estimate FODs
$(sd)/dwi/fod.mif: $(sd)/dwi/preproc.mif $(sd)/dwi/response.txt
	dwi2fod csd $^ $@ -force -nthreads $(nthread) $(dwi_log)

# convert FS labels to connectivity labels in DWI space
$(sd)/dwi/label.%.mif: $(sd)/dwi/aparc+aseg.%.nii.gz $(sd)/dwi/lut.%.txt
	labelconvert $< \
	    $(lut_fs) \
		$(sd)/dwi/lut.$*.txt \
	    $@ -force $(dwi_log)

# convert FS labels to connectivity labels in T1 space
.PRECIOUS: $(sd)/dwi/label_in_T1.%.nii.gz
$(sd)/dwi/label_in_T1.%.nii.gz: $(sd)/mri/aparc+aseg.%.RAS.RO.nii.gz $(sd)/dwi/lut.%.txt
	labelconvert $< \
		$(lut_fs) \
		$(sd)/dwi/lut.$*.txt \
		$@ -force $(dwi_log)

# generate all tracks
$(sd)/dwi/all.tck: $(sd)/dwi/fod.mif $(sd)/dwi/mask.mif
	tckgen $< $@ \
	    -mask $(sd)/dwi/mask.mif \
	    -seed_image $(sd)/dwi/mask.mif \
	    $(number_flag) $(ntrack) -force -nthreads $(nthread) $(dwi_log)

# subsample tracks
$(sd)/dwi/100k.tck: $(sd)/dwi/all.tck
	tckedit $(number_flag) 100K $< $@ $(dwi_log)

# generate track counts for connectome
$(sd)/dwi/triu_counts.%.txt: $(sd)/dwi/all.tck $(sd)/dwi/label.%.mif
	tck2connectome $^ $@ -force -nthreads $(nthread) $(dwi_log)

# generate track average lengths for connectome
$(sd)/dwi/triu_lengths.%.txt: $(sd)/dwi/all.tck $(sd)/dwi/label.%.mif
	tck2connectome $^ $@ -scale_length -stat_edge mean -force -nthreads $(nthread) $(dwi_log)

# convert to non-triangular, normalize, etc TODO
$(sd)/dwi/%.txt: $(sd)/dwi/triu_%.txt
	python -m util.util postprocess_connectome $< $@ $(dwi_log)


# Atlases
$(sd)/dwi/lut.dk.txt:
	cp `find $(MRT3) -name fs_default.txt | head -n 1` $@ $(dwi_log)

$(sd)/dwi/lut.destrieux.txt:
	cp `find $(MRT3) -name fs_a2009s.txt | head -n 1` $@ $(dwi_log)

$(sd)/dwi/lut.vep.txt:
	cp $(here)/util/data/VepMrtrixLut.txt $@ $(dwi_log)

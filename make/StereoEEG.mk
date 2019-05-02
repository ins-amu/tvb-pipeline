
# Coregister with T1
$(sd)/elec/elec-$(elec_mode)_in_T1.nii.gz: $(ELEC) $(fs_done) $(sd)/mri/T1.RAS.RO.nii.gz
	mkdir -p $(sd)/elec
	mri_convert $< $(sd)/elec/elec-$(elec_mode).nii.gz --out_orientation RAS
	flirt -in   $(sd)/elec/elec-$(elec_mode).nii.gz \
	      -omat $(sd)/elec/elec-$(elec_mode)_to_T1.mat \
	      -out  $(sd)/elec/elec-$(elec_mode)_in_T1.nii.gz \
          -ref  $(sd)/mri/T1.RAS.RO.nii.gz \
	    $(regopts)

$(sd)/elec/mask.nii.gz: $(sd)/mri/brain.RAS.nii.gz
	mkdir -p $(sd)/elec
	mri_binarize --i $< --o $@ --min 10 --erode 4

$(sd)/elec/masked_elec-CT.nii.gz: $(sd)/elec/elec-CT_in_T1.nii.gz $(sd)/elec/mask.nii.gz
	mri_binarize --i $< --o $@ --min $(ct_thresh) \
	    --mask $(sd)/elec/mask.nii.gz

$(sd)/elec/dilated_elec-CT.nii.gz: $(sd)/elec/masked_elec-CT.nii.gz
	mri_binarize --i $< --o $@ --min 0.5 --dilate 2 --erode 1

$(sd)/elec/labeled_elec-CT.nii.gz: $(sd)/elec/masked_elec-CT.nii.gz $(sd)/elec/dilated_elec-CT.nii.gz
	python -m util.util label_with_dilation $^ $@

# case of T1 w/ elec, no dilation
$(sd)/elec/masked_elec-T1.nii.gz: $(sd)/elec/elec-T1_in_T1.nii.gz $(sd)/elec/mask.nii.gz
	mri_binarize --i $< --o $@ --min $(elec_min) --max $(elec_max) --mask $(sd)/elec/mask.nii.gz

$(sd)/elec/labeled_elec-T1.nii.gz: $(sd)/elec/masked_elec-T1.nii.gz
	python -m util.util label_objects $< $@ $(nthread)


# wildcard syntax below tests if file exists
ifneq ("$(wildcard $(ELEC_POS_GARDEL))", "")
# Coordinates from GARDEL
$(sd)/elec/seeg.xyz: $(ELEC_POS_GARDEL) $(sd)/elec/elec-$(elec_mode)_in_T1.nii.gz
	python -m util.util transform_gardel_coords_to_tvb $< \
        $(sd)/elec/elec-$(elec_mode).nii.gz \
		$(sd)/elec/elec-$(elec_mode)_in_T1.nii.gz \
		$(sd)/elec/elec-$(elec_mode)_to_T1.mat $@
else ifneq ("$(wildcard $(ELEC_ENDPOINTS))", "")
# Marked endpoints of electrodes
$(sd)/elec/seeg.xyz: $(ELEC_ENDPOINTS) $(sd)/elec/elec-$(elec_mode)_in_T1.nii.gz
	python -m util.util gen_seeg_xyz_from_endpoints $< $@ \
	    $(sd)/elec/elec-$(elec_mode)_to_T1.mat \
		$(sd)/elec/elec-$(elec_mode).nii.gz	   \
		$(sd)/elec/elec-$(elec_mode)_in_T1.nii.gz
else
# User must label electrodes by hand at some point
$(sd)/elec/seeg.xyz: $(sd)/elec/labeled_elec-$(elec_mode).nii.gz $(ELEC_LABEL_SCHEMA)
	python -m util.util gen_seeg_xyz $^ $@
endif

$(sd)/elec/gain_dipole_no-subcort.%.txt: $(sd)/tvb/connectivity.%.zip $(sd)/elec/seeg.xyz
	python -m util.gain_matrix_seeg \
	  --mode surface \
	  --formula dipole \
	  --no_use_subcort \
	  --surf_dir $(sd)/tvb/ \
	  --parcellation $* \
	  $^ $@

$(sd)/elec/gain_inv-square.%.txt: $(sd)/tvb/connectivity.%.zip $(sd)/elec/seeg.xyz
	python -m util.gain_matrix_seeg \
	  --mode surface \
	  --formula inv_square \
	  --use_subcort \
	  --surf_dir $(sd)/tvb/ \
	  --parcellation $* \
	  $^ $@

$(sd)/elec/gain_volume.%.txt: $(sd)/dwi/label_in_T1.%.nii.gz $(sd)/tvb/connectivity.%.zip $(sd)/elec/seeg.xyz
	python -m util.gain_matrix_seeg             \
	  --mode volume                             \
	  --formula inv_square                      \
	  --use_subcort                             \
	  --label $(sd)/dwi/label_in_T1.$*.nii.gz   \
	  $(sd)/tvb/connectivity.$*.zip             \
      $(sd)/elec/seeg.xyz $@

$(sd)/elec/img: $(sd)/elec/seeg.xyz $(sd)/mri/T1.RAS.RO.nii.gz $(sd)/elec/elec-$(elec_mode)_in_T1.nii.gz
	mkdir -p $@
	python -m util.plot plot_t1_plus_elecs $(sd)/mri/T1.RAS.RO.nii.gz $(sd)/elec/elec-$(elec_mode)_in_T1.nii.gz \
	    $(sd)/elec/seeg.xyz $@

$(sd)/elec/elec.%.png: $(sd)/tvb/connectivity.%.zip $(sd)/elec/seeg.xyz
	python -m util.plot seeg_elecs $^ $@


# SEEG recordings -------------------------------------------------------------------- #
.PHONY: seeg
seeg: $(sd)/seeg/fif $(sd)/seeg/img

$(sd)/seeg/fif: $(XLSX) $(SEEGRECDIR) $(sd)/elec/seeg.xyz
	mkdir -p $@
	python -m util.parse_patient_xlsx convert_recordings $(XLSX) $(SEEGRECDIR) $(sd)/elec/seeg.xyz $(sd)/seeg/fif/

$(sd)/seeg/img: $(sd)/seeg/fif
	mkdir -p $@
	for i in `ls $(sd)/seeg/fif/*.json`; do \
		trg=$@/`basename $$i .json`; \
		python -m util.plot_seeg_recording $$i avgref $${trg}.avgref.png; \
		python -m util.plot_seeg_recording $$i bipolar $${trg}.bipolar.png; \
		python -m util.plot_seeg_recording $$i spectrogram $${trg}.spectrogram.png; \
	done;
# ------------------------------------------------------------------------------------ #

# Stereotactic EEG & CT

$(sd)/elec/CT_in_T1.nii.gz: $(CT) $(fs_done) $(sd)/mri/T1.RAS.nii.gz
	mkdir -p $(sd)/elec
	mri_convert $< $(sd)/elec/CT.nii.gz --out_orientation RAS
	flirt -in $(sd)/elec/CT.nii.gz -ref $(sd)/mri/T1.RAS.nii.gz \
	    -omat $(sd)/elec/CT_to_T1.mat \
	    -out $(sd)/elec/CT_in_T1.nii.gz \
	    $(regopts)

$(sd)/elec/mask.nii.gz: $(sd)/mri/brain.RAS.nii.gz
	mkdir -p $(sd)/elec
	mri_binarize --i $< --o $@ --min 10 --erode 4

$(sd)/elec/masked_CT.nii.gz: $(sd)/elec/CT_in_T1.nii.gz $(sd)/elec/mask.nii.gz
	mri_binarize --i $< --o $@ --min $(ct_thresh) \
	    --mask $(sd)/elec/mask.nii.gz

$(sd)/elec/dilated_CT.nii.gz: $(sd)/elec/masked_CT.nii.gz
	mri_binarize --i $< --o $@ --min 0.5 --dilate 2 --erode 1

$(sd)/elec/labeled_CT.nii.gz: $(sd)/elec/masked_CT.nii.gz $(sd)/elec/dilated_CT.nii.gz
	python -m util.util label_with_dilation $^ $@

# case of T1 w/ elec, no dilation
$(sd)/elec/ELEC_in_T1.nii.gz: $(ELEC) $(fs_done) $(sd)/mri/T1.RAS.nii.gz
	mkdir -p $(sd)/elec
	mri_convert $< $(sd)/elec/ELEC.nii.gz --out_orientation RAS
	flirt -in $(sd)/elec/ELEC.nii.gz -ref $(sd)/mri/T1.RAS.nii.gz \
	    -omat $(sd)/elec/ELEC_to_T1.mat \
	    -out $(sd)/elec/ELEC_in_T1.nii.gz \
	    $(regopts)

$(sd)/elec/masked_ELEC.nii.gz: $(sd)/elec/ELEC_in_T1.nii.gz $(sd)/elec/mask.nii.gz
	mri_binarize --i $< --o $@ --min $(elec_min) --max $(elec_max) --mask $(sd)/elec/mask.nii.gz

$(sd)/elec/labeled_ELEC.nii.gz: $(sd)/elec/masked_ELEC.nii.gz
	python -m util.util label_objects $< $@ $(nthread)

ifdef ELEC_ENDPOINTS
# more reliable to mark endpoints
$(sd)/elec/seeg.xyz: $(ELEC_ENDPOINTS) $(sd)/elec/$(elec_mode)_in_T1.nii.gz
	python -m util.util gen_seeg_xyz_from_endpoints $< $@ $(sd)/elec/$(elec_mode)_to_T1.mat \
	    $(sd)/elec/$(elec_mode).nii.gz $(sd)/elec/$(elec_mode)_in_T1.nii.gz
else ifdef ELEC_POS_GARDEL
$(sd)/elec/seeg.xyz: $(ELEC_POS_GARDEL) $(sd)/elec/$(elec_mode)_in_T1.nii.gz
	python -m util.util transform_gardel_coords_to_tvb $< \
        $(sd)/elec/$(elec_mode).nii.gz $(sd)/elec/$(elec_mode)_in_T1.nii.gz \
	    $(sd)/elec/$(elec_mode)_to_T1.mat $@
else
# user must label electrodes by hand at some point
$(sd)/elec/seeg.xyz: $(sd)/elec/labeled_$(elec_mode).nii.gz $(sd)/elec/schema.txt
	python -m util.util gen_seeg_xyz $^ $@
endif

$(sd)/elec/gain_dipole_no-subcort.mat: $(sd)/tvb/connectivity.zip $(sd)/elec/seeg.xyz
	python -m util.gain_matrix_seeg \
	  --mode surface \
	  --formula dipole \
	  --no_use_subcort \
	  --surf_dir $(sd)/tvb/ \
	  $^ $@

$(sd)/elec/gain_inv-square.mat: $(sd)/tvb/connectivity.zip $(sd)/elec/seeg.xyz
	python -m util.gain_matrix_seeg \
	  --mode surface \
	  --formula inv_square \
	  --use_subcort \
	  --surf_dir $(sd)/tvb/ \
	  $^ $@

$(sd)/elec/img: $(sd)/tvb/connectivity.zip $(sd)/elec/seeg.xyz $(sd)/mri/T1.RAS.nii.gz $(sd)/elec/$(elec_mode)_in_T1.nii.gz
	mkdir -p $@
	python -m util.plot seeg_elecs $(sd)/tvb/connectivity.zip $(sd)/elec/seeg.xyz $(sd)/elec/img/elec.png
	python -m util.plot plot_t1_plus_elecs $(sd)/mri/T1.RAS.nii.gz $(sd)/elec/$(elec_mode)_in_T1.nii.gz \
	    $(sd)/elec/seeg.xyz $@


# SEEG recordings -------------------------------------------------------------------- #
.PHONY: seeg
seeg: $(sd)/seeg/fif $(sd)/seeg/img

$(sd)/seeg/fif: $(XLSX) $(SEEGRECDIR)
	mkdir -p $@
	python -m util.parse_patient_xlsx convert_recordings $(XLSX) $(SEEGRECDIR) $(sd)/seeg/fif/

$(sd)/seeg/img: $(sd)/seeg/fif
	mkdir -p $@
	for i in `ls $(sd)/seeg/fif/*.json`; do \
		trg=$@/`basename $$i .json`; \
		python -m util.plot_seeg_recording $$i avgref $${trg}.avgref.png; \
		python -m util.plot_seeg_recording $$i bipolar $${trg}.bipolar.png; \
		python -m util.plot_seeg_recording $$i spectrogram $${trg}.spectrogram.png; \
	done;
# ------------------------------------------------------------------------------------ #

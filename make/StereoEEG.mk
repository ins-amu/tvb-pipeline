# Stereotactic EEG & CT

$(sd)/seeg/CT_in_T1.nii.gz: $(CT) $(fs_done) $(sd)/mri/T1.RAS.nii.gz
	mkdir -p $(sd)/seeg
	mri_convert $< $(sd)/seeg/CT.nii.gz --out_orientation RAS
	flirt -in $(sd)/seeg/CT.nii.gz -ref $(sd)/mri/T1.RAS.nii.gz \
	    -omat $(sd)/seeg/CT_to_T1.mat \
	    -out $(sd)/seeg/CT_in_T1.nii.gz \
	    $(regopts)

$(sd)/seeg/mask.nii.gz: $(sd)/mri/brain.RAS.nii.gz
	mkdir -p $(sd)/seeg
	mri_binarize --i $< --o $@ --min 10 --erode 4

$(sd)/seeg/masked_CT.nii.gz: $(sd)/seeg/CT_in_T1.nii.gz $(sd)/seeg/mask.nii.gz
	mri_binarize --i $< --o $@ --min $(ct_thresh) \
	    --mask $(sd)/seeg/mask.nii.gz

$(sd)/seeg/dilated_CT.nii.gz: $(sd)/seeg/masked_CT.nii.gz
	mri_binarize --i $< --o $@ --min 0.5 --dilate 2 --erode 1

$(sd)/seeg/labeled_CT.nii.gz: $(sd)/seeg/masked_CT.nii.gz $(sd)/seeg/dilated_CT.nii.gz
	python -m util.util label_with_dilation $^ $@

# case of T1 w/ elec, no dilation
$(sd)/seeg/ELEC_in_T1.nii.gz: $(ELEC) $(fs_done) $(sd)/mri/T1.RAS.nii.gz
	mkdir -p $(sd)/seeg
	mri_convert $< $(sd)/seeg/ELEC.nii.gz --out_orientation RAS
	flirt -in $(sd)/seeg/ELEC.nii.gz -ref $(sd)/mri/T1.RAS.nii.gz \
	    -omat $(sd)/seeg/ELEC_to_T1.mat \
	    -out $(sd)/seeg/ELEC_in_T1.nii.gz \
	    $(regopts)

$(sd)/seeg/masked_ELEC.nii.gz: $(sd)/seeg/ELEC_in_T1.nii.gz $(sd)/seeg/mask.nii.gz
	mri_binarize --i $< --o $@ --min $(elec_min) --max $(elec_max) --mask $(sd)/seeg/mask.nii.gz

$(sd)/seeg/labeled_ELEC.nii.gz: $(sd)/seeg/masked_ELEC.nii.gz
	python -m util.util label_objects $< $@ $(nthread)

ifdef ELEC_ENDPOINTS
# more reliable to mark endpoints
$(sd)/seeg/seeg.xyz: $(ELEC_ENDPOINTS) $(sd)/seeg/$(elec_mode)_in_T1.nii.gz
	python -m util.util gen_seeg_xyz_from_endpoints $< $@ $(sd)/seeg/$(elec_mode)_to_T1.mat \
	    $(sd)/seeg/$(elec_mode).nii.gz $(sd)/seeg/$(elec_mode)_in_T1.nii.gz
else
# user must label electrodes by hand at some point
$(sd)/seeg/seeg.xyz: $(sd)/seeg/labeled_$(elec_mode).nii.gz $(sd)/seeg/schema.txt
	python -m util.util gen_seeg_xyz $^ $@
endif

$(sd)/seeg/seeg.png: $(sd)/tvb/connectivity.zip $(sd)/seeg/seeg.xyz
	python -m util.plot seeg_elecs $^ $@

$(sd)/seeg/gain_dipole_no-subcort.mat: $(sd)/tvb/connectivity.zip $(sd)/seeg/seeg.xyz
	python -m util.gain_matrix_seeg \
	  --mode surface \
	  --formula dipole \
	  --no_use_subcort \
	  --surf_dir $(sd)/tvb/ \
	  $^ $@

$(sd)/seeg/gain_inv-square.mat: $(sd)/tvb/connectivity.zip $(sd)/seeg/seeg.xyz
	python -m util.gain_matrix_seeg \
	  --mode surface \
	  --formula inv_square \
	  --use_subcort \
	  --surf_dir $(sd)/tvb/ \
	  $^ $@


# SEEG recordings -------------------------------------------------------------------- #
.PHONY: seegrec
seegrec: $(sd)/seegrec/fif $(sd)/seegrec/img

$(sd)/seegrec/fif: $(XLSX) $(SEEGRECDIR)
	mkdir -p $@
	python $(here)/util/parse_patient_xlsx.py $(XLSX) $(sd)/seegrec/fif/
	for i in `find $(SEEGRECDIR) -maxdepth 1 -name '*.eeg'`; do \
		trg=$@/`basename "$${i%.*}"`.raw.fif; \
		python $(here)/util/read_eeg.py $$i $$trg; \
	done;

$(sd)/seegrec/img: $(sd)/seegrec/fif
	mkdir -p $@
	for i in `ls $(sd)/seegrec/fif/*.json`; do \
		trg=$@/`basename $$i .json`; \
		python $(here)/util/plot_seeg_recording.py $$i avgref $${trg}.avgref.png; \
		python $(here)/util/plot_seeg_recording.py $$i bipolar $${trg}.bipolar.png; \
		python $(here)/util/plot_seeg_recording.py $$i spectrogram $${trg}.spectrogram.png; \
	done;
# ------------------------------------------------------------------------------------ #

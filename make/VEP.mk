# prep data
$(sd)/vep/data.R: $(sd)/vep
	python -m util.vep_preprocess $(sd)

# ensure topic folder created once subject directory is available
$(sd)/vep: $(sd)/mri/orig/001.mgz
	mkdir -p $(sd)/vep

nchains := $(shell python -m util cores) 

# ensure topic folder created once subject directory is available
$(sd)/stan: $(sd)/mri/orig/001.mgz
	mkdir -p $(sd)/stan

# compile model
$(sd)/stan/%.stan.pkl: $(here)/stan/%.stan
	python -m util.stan compile_model $< $@

# generate data automatically if supporting py module present
$(sd)/stan/%.dat.pkl: $(here)/stan/%.dat.py
	python $< $@

# first optimize
$(sd)/stan/%.opt.pkl: $(sd)/stan/%.stan.pkl $(sd)/stan/%.dat.pkl
	python -m util.stan optimize_files $^ $@ $(nchains)

# then sample
$(sd)/stan/%.samp.pkl: $(sd)/stan/%.stan.pkl $(sd)/stan/%.dat.pkl \
                       $(sd)/stan/%.opt.pkl
	python -m util.stan sample_files $^ $@

# visualize?
$(sd)/stan/%.png: $(here)/stan/%.vis.py $(sd)/stan/%.samp.pkl
	python $^ $@

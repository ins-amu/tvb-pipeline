FROM ubuntu:20.04
LABEL maintainer="marmaduke.woodman@univ-amu.fr"

ENV DEBIAN_FRONTEND=noninteractive

WORKDIR /usr/local

RUN apt-get update \
 && apt-get install -y libnss3-dev libx11-xcb1 libxcb-dri3-0 libxcomposite1 \
	libxcursor1 libxdamage1 libxfixes3 libxi6 libxtst6 libatk1.0-0 libatk-bridge2.0-0 \
	libgdk-pixbuf2.0-0 libgtk-3-0 libgtk-3-0 libpangocairo-1.0-0 libpango-1.0-0 libcairo2 \
	libdrm2 libgbm1 libasound2 libatspi2.0-0 curl git build-essential tcsh perl nodejs \
	python2 wget datalad bc libglu1-mesa-dev unzip

RUN wget -q https://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/7.2.0/freesurfer-linux-ubuntu18_amd64-7.2.0.tar.gz \
 && tar xzf freesurfer-linux-ubuntu18_amd64-7.2.0.tar.gz \
 && rm freesurfer-linux-ubuntu18_amd64-7.2.0.tar.gz

RUN wget https://fsl.fmrib.ox.ac.uk/fsldownloads/fslinstaller.py \
 && echo "" | python2 fslinstaller.py

RUN curl -LO https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
 && bash Miniconda3-latest-Linux-x86_64.sh -b -p $PWD/conda \
 && rm Miniconda3-latest-Linux-x86_64.sh \
 && export PATH=$PWD/conda/bin:$PATH \
 && conda install -y jupyter numba scipy matplotlib \
 && conda install -y -c mrtrix3 mrtrix3 \
 && pip install tvb-data tvb-library tqdm pybids siibra requests pyunicore mne nilearn pyvista ipywidgets cmdstanpy \
 && install_cmdstan \
 && mv /root/.cmdstan $PWD/cmdstan

ENV PATH=/usr/local/conda/bin:$PATH
ENV FREESURFER_HOME=/usr/local/freesurfer

# TODO setup env automatically on entry

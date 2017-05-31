# To build this image, you'll need to download FreeSurfer, MNE, DCMTK
# and have a valid FreeSurfer license file.

FROM ubuntu:16.04
MAINTAINER Marmaduke Woodman <mmwoodman@gmail.com>

WORKDIR /opt

# system packages {{{
RUN apt-get update && apt-get install -y wget \
  && wget -O- http://neuro.debian.net/lists/xenial.de-m.full | tee /etc/apt/sources.list.d/neurodebian.sources.list \
  && apt-key adv --recv-keys --keyserver hkp://pool.sks-keyservers.net:80 0xA5D32F012649A5A9 \
  && apt-get update \
  && apt-get install -y g++ libeigen3-dev git python python-numpy zlib1g-dev \
	cmake tcsh libeigen3-dev liblapack-dev libblas-dev libssl-dev fsl-complete \
        libhdf5-dev zlib1g-dev libmatio-dev libopenblas-dev liblapacke-dev
# }}}

# external packages
ADD MNE*.tar.gz /opt/
ADD freesurfer*.tar.gz /opt/
ADD license.txt /opt/freesurfer/license.txt
ADD dcmtk*.tar.bz2 /opt
RUN mv dcmtk* dcmtk

# FS, FSL, MNE env vars {{{
ENV FIX_VERTEX_AREA= \
    FMRI_ANALYSIS_DIR=/opt/freesurfer/fsfast \
    FREESURFER_HOME=/opt/freesurfer \
    FSFAST_HOME=/opt/freesurfer/fsfast \
    FSF_OUTPUT_FORMAT=nii.gz \
    FS_OVERRIDE=0 \
    LOCAL_DIR=/opt/freesurfer/local \
    MINC_BIN_DIR=/opt/freesurfer/mni/bin \
    MINC_LIB_DIR=/opt/freesurfer/mni/lib \
    MNEROOT=/opt/MNE-2.7.0-3106-Linux-x86_64 \
    MNE_BIN_PATH=/opt/MNE-2.7.0-3106-Linux-x86_64/bin \
    MNE_LIB_PATH=/opt/MNE-2.7.0-3106-Linux-x86_64/lib \
    MNE_ROOT=/opt/MNE-2.7.0-3106-Linux-x86_64 \
    MNI_DATAPATH=/opt/freesurfer/mni/data \
    MNI_DIR=/opt/freesurfer/mni \
    MNI_PERL5LIB=/opt/freesurfer/mni/share/perl5 \
    OS=Linux \
    PERL5LIB=/opt/freesurfer/mni/share/perl5 \
    SUBJECTS_DIR=/opt/freesurfer/subjects \
    XUSERFILESEARCHPATH=/opt/MNE-2.7.0-3106-Linux-x86_64/share/app-defaults/%N \
    LD_LIBRARY_PATH=/opt/MNE-2.7.0-3106-Linux-x86_64/lib \
    PATH=/opt/MNE-2.7.0-3106-Linux-x86_64/bin:/opt/freesurfer/bin:/opt/freesurfer/fsfast/bin:/opt/freesurfer/tktools:/opt/freesurfer/mni/bin:/opt/dcmtk/bin:$PATH \
    DCMDICTPATH=/opt/dcmtk/share/dcmtk/dicom.dic \
    FSLDIR=/usr/share/fsl/5.0 \
    FSLBROWSER=/etc/alternatives/x-www-browser \
    FSLLOCKDIR= \
    FSLMACHINELIST= \
    FSLMULTIFILEQUIT=TRUE \
    FSLOUTPUTTYPE=NIFTI_GZ \
    FSLREMOTECALL= \
    FSLTCLSH=/usr/bin/tclsh \
    FSLWISH=/usr/bin/wish \
    LD_LIBRARY_PATH=/usr/lib/fsl/5.0:$LD_LIBRARY_PATH \
    PATH=/usr/lib/fsl/5.0:/usr/share/fsl/5.0/bin:$PATH \
    POSSUMDIR=/usr/share/fsl/5.0
# }}}

# Mrtrix3 {{{
RUN git clone https://github.com/mrtrix3/mrtrix3 \
 && cd mrtrix3 && ./configure -nogui && ./build
ENV PATH=/opt/mrtrix3/bin:$PATH
# }}}

# Openmeeg # {{{
RUN git clone https://github.com/openmeeg/openmeeg
RUN cd openmeeg/OpenMEEG && \
    mkdir build && cd build && \
    cmake -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release \
        -DENABLE_PYTHON=OFF -DCMAKE_INSTALL_PREFIX=/usr/local \
        -DBLASLAPACK_IMPLEMENTATION="OpenBLAS" \
        -DBUILD_DOCUMENTATION=OFF -DBUILD_TUTORIALS=OFF .. && \
    make -j && \
    make test && \
    make install
ENV LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
# }}}

# Py2 & TVB {{{
RUN curl -O https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
RUN bash Miniconda2-latest-Linux-x86_64.sh -b -p /opt/conda
ENV PATH /opt/conda/bin:$PATH
RUN conda install -c conda-forge jupyterlab \
 && conda install nomkl numpy numexpr numba matplotlib scipy cython scikit-learn \
 && pip install gdist psutil networkx nibabel nilearn mne
RUN for repo in library data; do git clone https://github.com/the-virtual-brain/tvb-$repo; done
ENV PYTHONPATH /opt/tvb-library:/opt/tvb-data:$PYTHONPATH
# }}}

# Ports {{{
EXPOSE 8888
EXPOSE 8000
# }}}

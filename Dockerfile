# To build this image, you'll need to download FreeSurfer, MNE
# and have a valid FreeSurfer license file, and put them into your scripts directory.
# To build:  cd path_to_scripts; docker build -t docker_scripts .
# To run: docker run -it -v <path_to_your_data_folder_on_host>:/opt/processing docker_scripts /bin/bash
# Then run scripts as usual: bash main_surface -c /opt/processing/<name_of_configuration_file.sh>

FROM ubuntu:16.04
MAINTAINER timpx <timpx@eml.cc>

# /opt used during installation, but 
# /opt/scripts is final workdir, set below
WORKDIR /opt

# system packages
RUN apt-get update && apt-get install -y wget \
  && wget -O- http://neuro.debian.net/lists/xenial.de-m.full | tee /etc/apt/sources.list.d/neurodebian.sources.list \
  && apt-get update \
  && apt-get install -y --allow-unauthenticated g++ libeigen3-dev git python \
      python-numpy zlib1g-dev fsl-complete gcc make cmake libopenblas-dev liblapacke-dev \
      libhdf5-dev libhdf5-serial-dev libmatio-dev python-dev swig libvtk6-dev doxygen libcgal-dev

# external packages
ADD freesurfer*.tar.gz /opt/
ADD license.txt /opt/freesurfer/license.txt

# FS, FSL, MNE env vars 
ENV FIX_VERTEX_AREA= \
    FREESURFER_HOME=/opt/freesurfer \
    FSFAST_HOME=/opt/freesurfer/fsfast \
    FSF_OUTPUT_FORMAT=nii.gz \
    FS_OVERRIDE=0 \
    LOCAL_DIR=/opt/freesurfer/local \
    MINC_BIN_DIR=/opt/freesurfer/mni/bin \
    MINC_LIB_DIR=/opt/freesurfer/mni/lib \
    SUBJECTS_DIR=/opt/freesurfer/subjects \
    MNI_DATAPATH=/opt/freesurfer/mni/data \
    MNI_DIR=/opt/freesurfer/mni \
    MNI_PERL5LIB=/opt/freesurfer/mni/share/perl5 \
    OS=Linux \
    PERL5LIB=/opt/freesurfer/mni/share/perl5 \
    MNE_ROOT=/opt/MNE-2.7.0-3106-Linux-x86_64 \
    MNE_BIN_PATH=/opt/MNE-2.7.0-3106-Linux-x86_64/bin \
    MNE_LIB_PATH=/opt/MNE-2.7.0-3106-Linux-x86_64/lib \
    MNE_ROOT=/opt/MNE-2.7.0-3106-Linux-x86_64 \
    XUSERFILESEARCHPATH=/opt/MNE-2.7.0-3106-Linux-x86_64/share/app-defaults/%N \
    LD_LIBRARY_PATH=/usr/lib/fsl/5.0:/opt/MNE-2.7.0-3106-Linux-x86_64/lib \
    PATH=/usr/lib/fsl/5.0:/usr/share/fsl/5.0/bin:/opt/MNE-2.7.0-3106-Linux-x86_64/bin:/opt/freesurfer/bin:/opt/freesurfer/fsfast/bin:/opt/freesurfer/tktools:/opt/freesurfer/mni/bin:/opt/dcmtk/bin:$PATH \
    FSLDIR=/usr/share/fsl/5.0 \
    FSLBROWSER=/etc/alternatives/x-www-browser \
    FSLLOCKDIR= \
    FSLMACHINELIST= \
    FSLMULTIFILEQUIT=TRUE \
    FSLOUTPUTTYPE=NIFTI_GZ \
    FSLREMOTECALL= \
    FSLTCLSH=/usr/bin/tclsh \
    FSLWISH=/usr/bin/wish \
    POSSUMDIR=/usr/share/fsl/5.0

# Mrtrix3
RUN git clone https://github.com/mrtrix3/mrtrix3 && cd mrtrix3 \
 && git checkout 0.3.16 && ./configure -nogui && ./build
ENV MRT3=/opt/mrtrix3
ENV PATH=/opt/mrtrix3/release/bin:/opt/mrtrix3/scripts:$PATH

# OpenMEEG
RUN git clone https://github.com/openmeeg/openmeeg \
 && cd openmeeg/OpenMEEG \
 && git checkout 2.4.prerelease \
 && mkdir build && cd build && \
    cmake -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release \
        -DENABLE_PYTHON=OFF -DCMAKE_INSTALL_PREFIX=/usr/local \
        -DBLASLAPACK_IMPLEMENTATION="OpenBLAS" \
        -DBUILD_DOCUMENTATION=OFF -DBUILD_TUTORIALS=OFF .. && \
    make -j && \
    make test && \
    make install
ENV LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH

# Py 3 stack
RUN curl -O https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda
ENV PATH=/opt/conda/bin:$PATH
RUN conda install -c conda-forge jupyterlab \
 && conda install numpy matplotlib \
 && pip install nibabel mne

# Makefile as entry point
RUN git clone https://github.com/ins-amu/scripts
WORKDIR /opt/scripts
RUN mkdir /opt/processing

#TODO for compatibility with tvb-make
# ENTRYPOINT ["make"] 





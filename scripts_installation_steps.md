The pipeline was tested on a Debian wheezy 64 bits, gcc version 4.6.3


    uname -a
    Linux prometheus 3.2.0-4-amd64 #1 SMP Debian 3.2.54-2 x86_64 GNU/Linux


Python
======

[Follow the snake:](http://continuum.io/downloads)

    which python
    /home/paula/anaconda/bin/python


Freesurfer
==========

[Download fressurfer Stable 5.3 for Linux x64:](http://surfer.nmr.mgh.harvard.edu/fswiki/DownloadAndInstall)
You will have to register and get a license key. 

[Install freesurfer](http://surfer.nmr.mgh.harvard.edu/fswiki/LinuxInstall)

    cd /usr/local/freesurfer

Add your license key:

    vim .license


Set the FREESURFER_HOME variable in .bashrc:

    export FREESURFER_HOME="/usr/local/freesurfer"


Configure freesurfer:

    paula@prometheus:~$ source $FREESURFER_HOME/SetUpFreeSurfer.sh
    -------- freesurfer-Linux-centos6_x86_64-stable-pub-v5.3.0 --------
    Setting up environment for FreeSurfer/FS-FAST (and FSL)
    FREESURFER_HOME   /usr/local/freesurfer
    FSFAST_HOME       /usr/local/freesurfer/fsfast
    FSF_OUTPUT_FORMAT nii.gz
    SUBJECTS_DIR      /usr/local/freesurfer/subjects
    INFO: /home/paula/matlab/startup.m does not exist ... creating
    MNI_DIR           /usr/local/freesurfer/mni


Change SUBJECTS_DIR and add to .bashrc
Create the subjects dir


Brainvisa
=========

[Download Brainvisa 4.4.0 for Linux 64 bits](http://brainvisa.info/downloadpage.html)

Install it:

    tar jxvf brainvisa-Mandriva-2008.0-x86_64-4.4.0-2013_11_18.tar.bz2 



Matlab
======

    Up to you

but you can install


Matlab Runtime Compiler
=======================

[Download MRC 8.1 for Matlab 2013a](wget http://www.mathworks.fr/supportfiles/MCR_Runtime/R2013a/MCR_R2013a_glnxa64_installer.zip)

    unzip MCR_R2013a_glnxa64_installer.zip
    ./intall -mode silent -agreeToLicense yes 


MRTRIX
======

[Download mrtrix 0.2-11](http://www.nitrc.org/projects/mrtrix)

Need an NTRIC account to download + license

[install mrtrix](http://www.brain.org.au/software/mrtrix/install/unix.html)

$ sudo apt-get install g++ python libglib2.0-dev libgtk2.0-dev libglibmm-2.4-dev libgtkmm-2.4-dev libgtkglext1-dev libgsl0-dev libgl1-mesa-dev libglu1-mesa-dev

                   tar xfj mrtrix-0.2.XX.bz2
                   cd mrtrix-0.2.12/
                   ./build 



Add to mrtrix to your PATH (in your .bashrc)

    /opt/mrtrix/bin

Before running the pipeline
    source $FREESURFER_HOME/SetUpFreeSurfer.sh

If you have installed mrtrix via `apt-get` from the Neurodebian repositories,
the binaries are installed to `/usr/lib/mrtrix/bin`.

FSL
===

You can download FSL on the [neurodebian repository](http://neuro.debian.net/pkgs/fsl-5.0-core.html?highlight=fsl)

to install neurodebian

    wget -O- http://neuro.debian.net/lists/wheezy.de-m.full | sudo tee /etc/apt/sources.list.d/neurodebian.sources.list
    sudo apt-key adv --recv-keys --keyserver pgp.mit.edu 2649A5A9

    sudo apt-get update
    
To install fsl-5.0-core:

    sudo apt-get install fsl-5.0-core










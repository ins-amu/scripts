# Pipeline 

#### Documentation
This pipeline uses T1 and diffusion imaging data to prepare surface, region mapping and connectivity for TVB in the following way:

- The surface is built with freesurfer
- The region mapping is built from the parcellation of freesurfer on a decimated surface (around 30.000 vertices)
- The connectivity and length matrices are build with mrtrix, i.e. spherical deconvolution and probabilistic tractrography. Right now there are two ways of building the connectivity matrix, either taking the ROIs only at the end of each fibres, or by taking into account all ROIs in which the fibres go through. We need further validation to say which one is the best one. The first method use around 100.000 tracks, the second around 400.000 tracks, however you can easily change these parameters.

#### Needed softwares
The pipeline was tested on a Debian wheezy 64 bits, gcc version 4.6.3
- freesurfer: tested with freesurfer-x86_64-redhat-linux-gnu-stable5-20130513
- brainvisa: tested with brainvisa 4.3.0
- python: tested with python 2.7.3
- matlab: tested with matlab R2013a
- mrtrix: tested with mrtrix 0.2.11
- fsl: tested with fsl 5.0

#### Run the pipeline 
- Create a main directory, all files will be in this main directory
In this main directory you must have:

 - scripts directory (clone of this repo)
 - data directory with:
    - T1 directory with T1 nifti or dicom
     - DWI directory with DWI data nifti or dicom

- Edit the four first lines (export) of scripts/main.sh to put the right directories instead.
Also edit the important parameters if needed.

- In a terminal:
```shell
cd main_directory/scripts
sh main.sh
```

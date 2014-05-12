# Surface and Connectivity Reconstruction: Imaging Pipeline & ToolS - SCRIPTS

#### Documentation
This pipeline uses T1 and diffusion imaging data to prepare surface, region mapping and connectivity for TVB in the following way:

- The surface is built with freesurfer
- The region mapping is built from the parcellation of freesurfer on a decimated surface (around 30.000 vertices)
- The connectivity and length matrices are build with mrtrix, i.e. spherical deconvolution and probabilistic tractrography. Right now there are two ways of building the connectivity matrix, either taking the ROIs only at the end of each fibres, or by taking into account all ROIs in which the fibres go through. We need further validation to say which one is the best one. The first method use around 100.000 tracks, the second around 400.000 tracks, however you can easily change these parameters.

Alternatively, you can run the region pipeline, which will provide only the necessary files for region simulation.
Conversely to the surface pipeline, you can use the parcellation you want with the region pipeline. To do so, put your parcellation nifti file and a correspondance matrix (which indicate the order of the region of your parcellation in the final weight matrix) in the parcellations directory.

#### Needed softwares
The pipeline was tested on a Debian wheezy 64 bits
- freesurfer: tested with freesurfer-x86_64-redhat-linux-gnu-stable5-20130513
- brainvisa: tested with brainvisa 4.3.0
- python: tested with python 2.7.3
- matlab or Matlab Compiler Runtime (free): tested with matlab R2013a and MCR 8.1
- mrtrix: tested with mrtrix 0.2.11
- fsl: tested with fsl 5.0

#### Installation
Please see [installation steps](scripts_installation_steps.md)

#### Run the pipeline 
- Create a main directory, all data and processed data files will be in this main directory
In this main directory you must have:

 - data directory with:
    - T1 directory with T1 nifti or dicom
    - DWI directory with DWI data nifti or dicom
    (note: if your DWI are in the nifti format, you will need to also provide the DW gradient scheme, with an extension in .b, such as explained in [Mrtrix documentation](http://www.brain.org.au/software/mrtrix/tractography/dwi.html))

- Clone this repository wherever you want, copy the example_config.sh file (you can copy it to your main directory) and edit it as needed.

- To run the surface pipeline, in a terminal:
```shell
cd path_to_scripts
bash main_surface.sh -c path_to_config/config.sh
```
To run the region pipeline, in a terminal:
```shell
cd path_to_scripts
bash main_region.sh -c path_to_config/config.sh
```

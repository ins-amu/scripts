# Surface and Connectivity Reconstruction: Imaging Pipeline for TVB Simulations
## SCRIPTS 0.4.1
 
Please see the [wiki](https://github.com/ins-amu/scripts/wiki)

#### New in 0.4.1:
- Remove need for Matlab
- new subparcellation scheme where arbitrary number of subregions can be chosen

#### New in 0.4:
- Dockerfile for easy installation with Docker.
- Tests. 
- Reorganization and simplification of the code => compatibility with tvb-make.
- 4 new parcellations: Destrieux, HCP-MMP1, Yeo7 and Yeo 17.

#### Features
- prepare data for TVB: surface reconstruction, region mapping, connectome.
- mrtrix 3.0 RC3: SIFT/SIFT2, ACT, 3 types of registration structural/diffusion, denoising, topup/eddy corrections, bias field corrections, mask upsampling and dilatation, multi-shell multi-tissue, dhollander algo, fsl subcortical structure parcellation, tractogram and tdi generation, multi-threaded.
- automatic config checks.
- convenience script for HCP datasets.
- handle automatically reverse phase-encoding DWI in most cases.
- subparcellation in any number of regions along with corresponding region mapping.
- MNE forward model.
- python 3.5.

#### License
This poject use the MIT License.
The full license is in LICENSE.txt in the SCRIPTS distribution.

Copyright (c) [2014] [The SCRIPTS Developers]

#### How to cite

When citing SCRIPTS, please cite this work:

Proix T, Spiegler A, Schirner M, Rothmeier S, Ritter P, Jirsa, VK. How do parcellation size and short-range connectivity affect dynamics in large-scale brain network models? NeuroImage, 2016, 142:135-149

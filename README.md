# Surface and Connectivity Reconstruction: Imaging Pipeline for TVB Simulations
## SCRIPTS 0.4-CH
 
Please see the [wiki](https://github.com/ins-amu/scripts/wiki)

If using SCRIPTS, please cite this work:
Proix T, Spiegler A, Schirner M, Rothmeier S, Ritter P, Jirsa, VK. How do parcellation size and short-range connectivity affect dynamics in large-scale brain network models? NeuroImage, 2016, 142:135-149 \[[link](https://www.sciencedirect.com/science/article/abs/pii/S1053811916302518)\]

If using the Connectomer Harmonics add-on, please cite:
Naze S., Proix T., Atasoy S., Kozloski J. R. (2020) Robustness of connectome harmonics to local gray matter and long-range white matter connectivity changes. NeuroImage \[[link](https://www.sciencedirect.com/science/article/pii/S1053811920308508)\]


#### Connectome Harmonics release:
- Compute high-resolution connectome based on tracks - surface mesh intersections (requires MATLAB Parallel Computing Toolbox).
- Compute eigenvalue/eigenvector pairs of combined local and long-range connectivity matrices.
- Compute Mutual Information of 100 harmonics (of lowest eigenvalues) with the Default Mode Network.
- To perform: after running `main_surface.sh`, run `main_CH.sh`.
- Template connectome harmonics for surfaces *cvs\_avg35\_inMNI152*, *fsaverage4* and *fsaverage5* using the [Gibbs tractography streamlines](https://www.nitrc.org/projects/gibbsconnectome) is now [available on zenodo](https://zenodo.org/record/4027989)


#### Connectome Harmonics dependencies
Install the following packages in the matlab_utils folder:
- Toolbox graph <https://www.mathworks.com/matlabcentral/fileexchange/5355-toolbox-graph>
- CBrewer colormaps  <https://www.mathworks.com/matlabcentral/fileexchange/34087-cbrewer---colorbrewer-schemes-for-matlab>
- Export_fig <https://www.mathworks.com/matlabcentral/fileexchange/23629-export-fig>
- TriangleRayIntersection <https://www.mathworks.com/matlabcentral/fileexchange/33073-triangle-ray-intersection>
- SmoothPatch <https://www.mathworks.com/matlabcentral/fileexchange/26710-smooth-triangulated-mesh>
- Subtightplot <https://www.mathworks.com/matlabcentral/fileexchange/39664-subtightplot>

#### NEW in 0.4:
- Dockerfile for easy installation with Docker.
- Tests. 
- Reorganization and simplification of the code => compatibility with tvb-make.
- 2 new parcellations: Destrieux, and HCP-MMP1

#### Features
- prepare data for TVB: surface reconstruction, region mapping, connectome.
- mrtrix 3.0 RC2: SIFT/SIFT2, ACT, 3 types of registration structural/diffusion, denoising, topup/eddy corrections, bias field corrections, mask upsampling and dilatation, multi-shell multi-tissue, dhollander algo, fsl subcortical structure parcellation, tractogram and tdi generation, multi-threaded.
- automatic config checks.
- convenience script for HCP datasets.
- handle automatically reverse phase-encoding DWI in most cases.
- subparcellation in any number of regions along with corresponding region mapping.
- MNE forward model.
- python 3.5.

#### License
This poject use the MIT License.
The full license is in LICENSE.txt in the SCRIPTS distribution.

Copyright (c) [2014-2020] [The SCRIPTS Developers]



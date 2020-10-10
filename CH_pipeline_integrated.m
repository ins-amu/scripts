% Script to compute connectome harmonics
% 	Processing steps: 
%		1) Import cortical surface surface mesh
%		2) Import tracks
%		3) Compute high-resolution connectome
%		4) Compute harmonics
%		5) Compute mutual information with default mode network 
%
% Author: Sebastien Naze
%
% Dependencies (none are included in repo, those with + sign must be added in matlab_utils/): 
%   - Gibbs Connectomne (https://www.nitrc.org/frs/?group_id=987)
%   - Freesurfer cortical surface templates (https://surfer.nmr.mgh.harvard.edu/fswiki/DownloadAndInstall)
%   + and surface readers 
%   + Toolbox graph (https://www.mathworks.com/matlabcentral/fileexchange/5355-toolbox-graph)
%   + CBrewer colormaps  (https://www.mathworks.com/matlabcentral/fileexchange/34087-cbrewer---colorbrewer-schemes-for-matlab)
%   + Export_fig (https://www.mathworks.com/matlabcentral/fileexchange/23629-export-fig)
%   + TriangleRayIntersection (https://www.mathworks.com/matlabcentral/fileexchange/33073-triangle-ray-intersection)
%   + SmoothPatch (https://www.mathworks.com/matlabcentral/fileexchange/26710-smooth-triangulated-mesh)
%   + Subtightplot (https://www.mathworks.com/matlabcentral/fileexchange/39664-subtightplot)

% global parameters
global FS PRD SUBJ_ID
FS = getenv('SUBJECTS_DIR');
PRD = getenv('PRD');
SUBJ_ID = getenv('SUBJ_ID');
nWorkers = str2num(getenv('NWORKERS'));
SCRIPTS_DIR = getenv('SCRIPTS_DIR');

% path to dependencies
addpath(genpath(fullfile(SCRIPTS_DIR, 'matlab_utils')));

% Flags for cortical surface template mesh to use (if using template)
CVS_Avg35_inMNI152 = 0;
fsaverage4 = 0;
fsaverage5 = 0;
fsaverage6 = 0;

% Flag for subject-specific cortical surface
indiv_surf = 1;

% Flag for using pial or white matter surface (default: WM)
pial = 0; 

% Flags for cortical surface mesh smoothing
smoothIter = 0;  % <-- keep at 0, not available in public realease. Contact author.  

% Flags for tractography
Yeh2018_dsistudio = 0;
GibbsConnectome = 0;
indiv_tracks = 1;
reload_tracks = 1;

% Flags for connectome
compute_connectome = 1;
save_connectome = 1;
plot_connectome = 0;
plot_tckLengths = 0;
compute_tck_stats = 0;
plot_tck_stats = 0;
plot_connectomeStats = 0; % <-- needs add-on, contact author 

% Parameters for connectome
nsamples = 0;                   % 0: all tracks, otherwise nbr of tracks
sz_bound = 3;			% parameter k_s in paper
ext_scaling = 2; 		% length of the track extension in k_s unit
minTckLength = 10; 		% ~ 1mm per unit, depending on tractography algo

% Flags and parameters for connectome harmonics
compute_harmonics = 1;
reload_connectome = 0;
nHarmonics = 100;		% Number of harmonics to compute
threshold_adjacency = -99;	% z_C in paper (use negative value << 0 to include all tracks)
trim_adjacency = 0;  % trimming
rand_trim_type = 'Ascend';
tckLen_trim_adjacency = 0;  
tckLen_trim_type = '';
local_spread = 1;		% nbr of neighbooring vertices for local diffusion on mesh
anisotropy = 0;			% percentage of local connectivity trimming (max 1)
  dist_based_local_trim = 0;  % distance based anisotropy, set to 1 to perform
    AniSorting = 'Descend';   % remove connections by ascending or descending order of length
callosectomy = 0;   % percentage of inter-hemispheric connections trimmed (max 1)
  n_graph_chunks = 1;   % number of subgraphs to perform the eigendecomposition on (set to 2 if disconnected hemispheres)
  dist_based_callo_trim = 1;  % distance-based callosectomy
    CalloSorting = 'Ascend';  % remove connections by ascending or descending order of length
separate_hemispheres = 0;     % set to 1 if disconnected hemispheres so that both are shown when projected on cortical surface
plot_CC_stats = 0;
rand_graph = [0 0 0];		% Flag for randomization type [interHemis, intraHemis, global]]
randomize_graph_CC_only = rand_graph(1);
randomize_graph_intraHemi_only = rand_graph(2);
randomize_graph_global = rand_graph(3);
rand_proba = 1;   % for performing gradual randomization (default:1)
iter_rand = 0;    % for numbering surrogates
sigma_eigs = 10^-9; % param of eigendec, check help eig
tol_eigs = 10^-9;   %   "
plot_graph = 0;     % set to 1 to plot graph illustration
plot_harmonics = 0; % set to 1 to show harmonics projected one by one on cortical surface (need plotNFsurf add-on, contact author)
nPlottedHarmonics = 20;		% nbr of harmonics to visualize

% Flags and params for MI DMN
MI_DMN = 1;
RSN='DMN';
RSN_Yeo2011 = 0;
reload_harmonics = 0;
nBinsMI = [16]; % always keep 1 first (if v3 wanted) and ascending order
thZscMI = 1;
MI_ylim = [0 0.3];
save_MI = 1;

% Flags and params for MI DSK
convert_into_DSK_atlas = 1; % /!\ (make sure to reload harmonics if want to load combined)
plot_converted_harmonics = 0;
re_save_combined = 1;

% Global flags for visualization and savings
vis = 'off';
plot_DMN_mask = 0; % <-- needs plotNFsurf add-on
save_combined = 1;
save_figs = 1;
close_figs = 1;

% miscallenous parameters
seed_rng = 'shuffle'; %12345; % replace by 'shuffle' to get seed based on current time
rng(seed_rng);
%nWorkers = 4;      % nbr of process for parallel computations (needs parallel processing toolbox) - commented because given as environment variable
nz = ceil(nsamples/nWorkers/2);
subj_type = 'HCP';	% subject type, only used for colormaps. HCP: Blue-White-Red; HD: Purple-White-Green

if randomize_graph_global
    rand_type = strcat('_randGlobal', num2str(iter_rand));
elseif randomize_graph_CC_only
    rand_type = strcat('_randCC', num2str(iter_rand));
    if randomize_graph_intraHemi_only
        rand_type = strcat('_randIntraAndCC', num2str(iter_rand));
    end
elseif randomize_graph_intraHemi_only
    rand_type = strcat('_randIntra', num2str(iter_rand));
else    
    if iter_rand>0
      rand_type = strcat('_rand', num2str(iter_rand));
    else
      rand_type = '';
    end
end

if rand_proba==1
  rnd_proba = '';
else
  rnd_proba = strcat('_rndP', num2str(rand_proba*100));
end

system(strcat("mkdir -p ", PRD, '/connectivity/img_', SUBJ_ID));

if trim_adjacency > 0
    trimAdj_prefix = strcat('_randTrim', rand_trim_type);
    trimAdj_all = strcat(trimAdj_prefix, num2str(trim_adjacency*100));
elseif tckLen_trim_adjacency > 0
    trimAdj_prefix = strcat('_tckLenTrim', tckLen_trim_type);
    trimAdj_all = strcat(trimAdj_prefix, num2str(tckLen_trim_adjacency*100));
else 
    trimAdj_prefix = '';
    trimAdj_all = '';
end

if RSN_Yeo2011
  rsn_type = strcat(RSN, '_Yeo2011_');
else
  rsn_type = strcat(RSN, '_rsn_');
end

if n_graph_chunks > 1
  callo_suffix = strcat('_nchunks', num2str(n_graph_chunks));
else
  callo_suffix = '';
end

if callosectomy == 1
  %separate_hemispheres = 1;
end

if separate_hemispheres
  sep_hemi='_sepHemiCombi';
else
  sep_hemi='';
end

if dist_based_local_trim
  distAni = strcat('DistBased', AniSorting(1:3));
else
  distAni = '';
end

if dist_based_callo_trim
  distCallo = CalloSorting(1:3);
else
  distCallo = '';
end

%% IMPORT SURFACE %%
if CVS_Avg35_inMNI152
    SUBJ_ID = 'copy_cvs_avg35_inMNI152';
    PRD = strcat(FS, SUBJ_ID);
    surf_type = 'cvsAvgInMNI152';
    
    if pial
        surf_type = strcat(surf_type, '_pial');
        vertices.left = load(strcat(PRD,'/surface/lh_vertices_low.txt'));
        vertices.right = load(strcat(PRD,'/surface/rh_vertices_low.txt'));
        
        faces.left = load(strcat(PRD,'/surface/lh_triangles_low.txt'))+1;  % +1 because matlab indexing starts at 1
        faces.right = load(strcat(PRD,'/surface/rh_triangles_low.txt'))+1;  % +1 because matlab indexing starts at 1
        
        vertices.all = [vertices.left; vertices.right];
        faces.all = [faces.left; faces.right + size(vertices.left,1)];
        
        % read white matter region mapping
        wm_regmap_lh = load(strcat(PRD, '/surface/lh_region_mapping_low.txt'));
        wm_regmap_rh = load(strcat(PRD, '/surface/rh_region_mapping_low.txt'));
        wm_regmap = [wm_regmap_lh; wm_regmap_rh];
        
        % considers smoothed surface only from white
    elseif ~smoothIter
        % read white matter surface
        vertices.left = load(strcat(PRD,'/surface/lh_white_vertices_low.txt'));
        vertices.right = load(strcat(PRD,'/surface/rh_white_vertices_low.txt'));
        faces.left = load(strcat(PRD,'/surface/lh_white_triangles_low.txt'))+1;  % +1 because matlab indexing starts at 1
        faces.right = load(strcat(PRD,'/surface/rh_white_triangles_low.txt'))+1;  % +1 because matlab indexing starts at 1
        vertices.all = [vertices.left; vertices.right];
        faces.all = [faces.left; faces.right + size(vertices.left,1)];
        
        % read white matter region mapping
        wm_regmap_lh = load(strcat(PRD, '/surface/lh_white_region_mapping_low.txt'));
        wm_regmap_rh = load(strcat(PRD, '/surface/rh_white_region_mapping_low.txt'));
        wm_regmap = [wm_regmap_lh; wm_regmap_rh];
    else
        % import surface
        surf_path = fullfile(PRD, 'surface', strcat('smoothSurf_iter', num2str(smoothIter)));
        vertices.left = dlmread(strcat(surf_path, '_lh_vertices.txt'), ' ');
        vertices.right = dlmread(strcat(surf_path, '_rh_vertices.txt'), ' ');
        %vertices.all = [vertices.left; vertices.right];
        faces.left = dlmread(strcat(surf_path, '_lh_faces.txt'), ' ');
        faces.right = dlmread(strcat(surf_path, '_rh_faces.txt'), ' ');
        vertices.all = [vertices.left; vertices.right];
        faces.all = [faces.left; faces.right + size(vertices.left,1)];
        
    end
    
elseif (fsaverage4 || fsaverage5 || fsaverage6)
    if fsaverage4
        SUBJ_ID ='copy_fsaverage4';
        surf_type = 'fsaverage4';
        disp('Importing fsaverage4 cortical surface...');
    elseif fsaverage5
        SUBJ_ID ='copy_fsaverage5';
        surf_type = 'fsaverage5';
        disp('Importing fsaverage5 cortical surface...');
    elseif fsaverage6
        SUBJ_ID ='copy_fsaverage6';
        surf_type = 'fsaverage6';
        disp('Importing fsaverage5 cortical surface...');
    end
    PRD = strcat(FS, SUBJ_ID);
    
    if pial
        [vertices.left, faces.left] = freesurfer_read_surf(strcat(PRD, '/surf/lh.pial'));
        [vertices.right, faces.right] = freesurfer_read_surf(strcat(PRD, '/surf/rh.pial'));   
        surf_type = strcat(surf_type, '_pial');
    else
        [vertices.left, faces.left] = freesurfer_read_surf(strcat(PRD, '/surf/lh.white'));
        [vertices.right, faces.right] = freesurfer_read_surf(strcat(PRD, '/surf/rh.white'));
        surf_type = strcat(surf_type, '_white');
    end
    
    vertices.all = [vertices.left; vertices.right];
    faces.all = [faces.left; faces.right + size(vertices.left,1)];
    
    [v.left, l.left, c.left] = read_annotation(strcat(FS, SUBJ_ID,'/label/lh.aparc.annot'));
    [v.right, l.right, c.right] = read_annotation(strcat(FS, SUBJ_ID,'/label/rh.aparc.annot'));
    lhreftable = load(strcat(SCRIPTS_DIR, 'util/lh_ref_table.txt', '-ascii'));
    rhreftable = load(strcat(SCRIPTS_DIR, 'util/rh_ref_table.txt', '-ascii'));
    
    % convert region mapping from fs to scripts
    wm_regmap_lh = l.left;
    for i = 1:numel(wm_regmap_lh)
        ref_idx = find(lhreftable(:,6) == wm_regmap_lh(i));
        if ~isempty(lhreftable(ref_idx,5))
            wm_regmap_lh(i) = lhreftable(ref_idx,5);
        else 
            wm_regmap_lh(i) = 0;
        end
    end
    wm_regmap_rh = l.right;
    for i = 1:numel(wm_regmap_rh)
        ref_idx = find(rhreftable(:,6) == wm_regmap_rh(i));
        if ~isempty(lhreftable(ref_idx,5))
            wm_regmap_rh(i) = rhreftable(ref_idx,5);
        else
            wm_regmap_rh(i) = 0;
        end
    end
     wm_regmap = [wm_regmap_lh; wm_regmap_rh];

     % Correct for temporal poles
     temp = zeros(size(vertices.all,1),1);
     unknown_idx = find(wm_regmap==0);
     tpoles_idx = unknown_idx(vertices.all(unknown_idx,3)<-15);
     wm_regmap(tpoles_idx) = 0;
     lh_tpoles_idx = tpoles_idx(tpoles_idx <= size(vertices.left,1));
     rh_tpoles_idx = tpoles_idx(tpoles_idx > size(vertices.left,1));
     wm_regmap(lh_tpoles_idx) = 42;
     wm_regmap(rh_tpoles_idx) = 85;
    
    
elseif indiv_surf
    surf_type = SUBJ_ID;
    if pial
        vertices.left = dlmread(strcat(PRD, '/surface/lh_vertices_low.txt'), ' ');
        vertices.right = dlmread(strcat(PRD, '/surface/rh_vertices_low.txt'), ' ');
        vertices.all = [vertices.left; vertices.right];
        faces.left = load(strcat(PRD,'/surface/lh_triangles_low.txt'))+1;  % +1 because matlab indexing starts at 1
        faces.right = load(strcat(PRD,'/surface/rh_triangles_low.txt'))+1;  % +1 because matlab indexing starts at 1
        faces.all = [faces.left; faces.right + size(vertices.left,1)];

        % read white matter region mapping
        wm_regmap_lh = load(strcat(PRD, '/surface/lh_region_mapping_low.txt'));
        wm_regmap_rh = load(strcat(PRD, '/surface/rh_region_mapping_low.txt'));
        wm_regmap = [wm_regmap_lh; wm_regmap_rh];
    else
        vertices.left = dlmread(strcat(PRD, '/surface/lh_white_vertices_low.txt'), ' ');
        vertices.right = dlmread(strcat(PRD, '/surface/rh_white_vertices_low.txt'), ' ');
        vertices.all = [vertices.left; vertices.right];
        faces.left = load(strcat(PRD,'/surface/lh_white_triangles_low.txt'))+1;  % +1 because matlab indexing starts at 1
        faces.right = load(strcat(PRD,'/surface/rh_white_triangles_low.txt'))+1;  % +1 because matlab indexing starts at 1
        faces.all = [faces.left; faces.right + size(vertices.left,1)];

        % read white matter region mapping
        wm_regmap_lh = load(strcat(PRD, '/surface/lh_white_region_mapping_low.txt'));
        wm_regmap_rh = load(strcat(PRD, '/surface/rh_white_region_mapping_low.txt'));
        wm_regmap = [wm_regmap_lh; wm_regmap_rh];
    end
end


%%%%%%%%%%%%%%%%%%%
%% IMPORT TRACKS %%
%%%%%%%%%%%%%%%%%%%
if GibbsConnectome
    connectome_type = 'Gibbs';
    % import tracks
    if reload_tracks
        load('/your/path/to/groupconnectome_169_horn2015.mat');	% to download, see HEADER
        ntcks = numel(fibs);
    end
    
elseif indiv_tracks
    connectome_type = SUBJ_ID;
    
    brain_info = niftiinfo(strcat(PRD, '/connectivity/brain.nii.gz'));
    M_vox2mm = brain_info.Transform.T(1:3, 1:3);
    A_vox2mm = brain_info.Transform.T(end, 1:3);
    
    M_RAS2LAS = [[-1,0,0]; [0,1,0]; [0,0,1]];
    A_RAS2LAS = [brain_info.ImageSize(1), 0, 0];
    
    system(strcat("mkdir -p ", PRD, '/connectivity/tmp_ascii_tck'));

    if (reload_tracks || (compute_connectome && ~exist('fibs', 'var')))
        files = dir(strcat(PRD, '/connectivity/tmp_ascii_tck/*.txt'));
        ntcks = numel(files);
        tck_idx = 1:ntcks;
        tck_idx = tck_idx(randperm(numel(tck_idx)));
        fibs = {};
        if ~nsamples
            nsamples = ntcks;
        else
            ntcks = nsamples;
        end
        tckidx = 1;
        for i = 1:nsamples
            tck_i = load(strcat(PRD, '/connectivity/tmp_ascii_tck/',files(i).name));
            if size(tck_i, 1) >= minTckLength
                fibs{tckidx} = (((tck_i * M_RAS2LAS) + A_RAS2LAS) * M_vox2mm) + A_vox2mm;
                tckidx = tckidx + 1;
            end
            if ~mod(i,10000)
                fprintf('%i percent of tracks imported \n\r', round(i*100/nsamples));
            end
        end
        ntcks = tckidx-1;
        if nsamples > ntcks
            nsamples = ntcks;
        end
    else
        ntcks = nsamples
    end
    
end

     
% update suffix for saving files
save_suffix = strcat('_', surf_type, '_', connectome_type, '_nsamples', num2str(nsamples), '_smoothIter', num2str(smoothIter), '_szBnd', num2str(sz_bound), '_extScaling', num2str(ext_scaling), '_minTckLen', num2str(minTckLength)); 

%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute connectome %%
%%%%%%%%%%%%%%%%%%%%%%%%
if compute_connectome
    
    % tracks stats    
    tck_idx = 1:ntcks;
    tck_idx = tck_idx(randperm(numel(tck_idx)));
    
    % tracks to parse
    if nsamples
        tck_idx = tck_idx(1:nsamples);
    end
    
    % setting up parallel pools
    nTcksPerWorker = floor(nsamples/nWorkers);
    workersIdx = 0:nTcksPerWorker:nsamples;
    parpool(nWorkers);
    dfibs = distributed(fibs(tck_idx));
    spmd
        % allocate memory of individual workers
        localConnectome = sparse([], [], [], size(vertices.all, 1), size(vertices.all, 1), nz);
        local_tckLengths_nb = sparse([], [], [], size(vertices.all, 1), size(vertices.all, 1), nz);
        
        % tracks stats
        if compute_tck_stats
            n_isects_w = [];  % number of intersections with cortical surface mesh
            stepSz_w = [];    % distance between 2 consecutive pts
            tckLengths_mm_w = [];
        end
    
        own_fibs = getLocalPart(dfibs);
        for i = 1:numel(own_fibs) 
            tck = own_fibs{i};

            % display number of tracks processed
            if ~mod(i,1000)
                fprintf('%i tcks done out of %i for worker %i \n', i, nTcksPerWorker, labindex);
            end

            % if track long enough for processing, proceed
            if ( (size(tck,1) > sz_bound*2) && (size(tck,1) > minTckLength) ) 
                % first bound
                bound1 = tck(1:sz_bound,:);
                b1vec = ext_scaling * (bound1(end,:) - bound1(1,:));
                [isect1, ~, ~, ~, xcoords1] = TriangleRayIntersection(bound1(end,:), bound1(end,:) + b1vec, vertices.all(faces.all(:,1),:), vertices.all(faces.all(:,2),:), vertices.all(faces.all(:,3),:), 'lineType', 'segment');

                if compute_tck_stats
                    n_isects_w = [sum(isect1), n_isects_w]; 
                end

                % check second bound only if first bound is valid
                if sum(isect1) > 0
                    % if several crossings, take the one closest to the beginning of the bound segment
                    if sum(isect1) > 1
                        dists1 = pdist2(xcoords1, bound1(end,:));
                        [minVal minIdx] = min(dists1);
                        xcoord1 = xcoords1(minIdx,:);    
                        isect1_face_idx = minIdx;
                    else
                        isect1_face_idx = find(isect1);
                        xcoord1 = xcoords1(isect1_face_idx, :);
                    end

                    % second bound
                    bound2 = tck(end-sz_bound:end,:);
                    b2vec = ext_scaling * (bound2(1,:) - bound2(end,:));
                    [isect2, ~, ~, ~, xcoords2] = TriangleRayIntersection(bound2(1,:), bound2(1,:) + b2vec, vertices.all(faces.all(:,1),:), vertices.all(faces.all(:,2),:), vertices.all(faces.all(:,3),:), 'lineType', 'segment');

                    if compute_tck_stats
                        n_isects_w = [sum(isect2), n_isects_w]; 
                    end

                    if sum(isect2) > 0
                        % if several crossings, take the one closest to the beginning of the bound segment
                        if sum(isect2) > 1
                            dists2 = pdist2(xcoords2, bound2(1,:));
                            [minVal minIdx] = min(dists2);
                            xcoord2 = xcoords2(minIdx,:);    
                            isect2_face_idx = minIdx;
                        else
                            isect2_face_idx = find(isect2);
                            xcoord2 = xcoords2(isect2_face_idx, :);
                        end

                        % get for each xcoords the corresponding closest vertex from the intersecting faces
                        dists1 = pdist2(xcoord1, vertices.all(faces.all(isect1_face_idx,:),:));
                        [minDist1 minIdx1] = min(dists1);
                        dists2 = pdist2(xcoord2, vertices.all(faces.all(isect2_face_idx,:),:));
                        [minDist2 minIdx2] = min(dists2);

                        A = faces.all(isect1_face_idx,minIdx1);
                        B = faces.all(isect2_face_idx,minIdx2);

                        % add connection to connectome
                        localConnectome(A,B) = localConnectome(A, B) + 1;
    
                        % add nbr of step to tckLengths..
                        local_tckLengths_nb(A,B) = local_tckLengths_nb(A,B) + size(tck,1);

                        % add step size and track length to stats
                        if compute_tck_stats
                            d = diff(tck);
                            steps = sqrt(sum(d.*d,2));
                            stepSz_w = [stepSz_w, steps'];
                            tckLengths_mm_w = [sum(steps), tckLengths_mm_w];
                        end
                    end
                end
            end
        end
    end
    
    % global memory allocation
    connectome = sparse(size(vertices.all, 1), size(vertices.all, 1));
    tckLengths_nb = sparse(size(vertices.all, 1), size(vertices.all, 1));
    n_isects = [];  % number of intersections with cortical surface mesh
    stepSz = [];    % distance between 2 consecutive pts
    tckLengths_mm = [];
    
    for ii = 1:nWorkers
        connectome = connectome + localConnectome{ii};
        tckLengths_nb = tckLengths_nb + local_tckLengths_nb{ii};
        if compute_tck_stats
            n_isects = [n_isects, n_isects_w{ii}];
            stepSz = [stepSz, stepSz_w{ii}];
            tckLengths_mm = [tckLengths_mm, tckLengths_mm_w{ii}];
        end
    end    
    delete(gcp);
    
    tckLengths_nb = (tckLengths_nb + tckLengths_nb') ./ (connectome + connectome');
    tckLengths_nb(isnan(tckLengths_nb)) = 0;
    
    if save_connectome
        save(strcat(PRD, '/connectivity/highRes_connectome', save_suffix, '.mat'), 'connectome');
        save(strcat(PRD, '/connectivity/highRes_tckSize', save_suffix, '.mat'), 'tckLengths_nb');
    end
    
    if plot_connectome
        connectome_fig = figure('Position', [600 600 1200 600], 'Visible', vis);
        subplot(1,2,1);
        spy(connectome + connectome');
        title('Connectome');
        
        subplot(1,2,2);
        imagesc(tckLengths_nb);
        title('Track lengths (nbr of segments)');
        
        if save_figs
            fname = strcat('highRes_connectivityMatrix_', connectome_type, '_smoothIt', num2str(smoothIter), '_nsamples', num2str(nsamples), '_szBnd', num2str(sz_bound), '_extScaling', num2str(ext_scaling));
            pause(2)
            export_fig(strcat(PRD, '/connectivity/img_',SUBJ_ID,'/', fname), '-png', '-transparent');
            pause(2);
            export_fig(strcat(PRD, '/connectivity/img_',SUBJ_ID,'/', fname), '-eps');
            pause(2)
            display(strcat(fname, ' figure exported'));
            if close_figs
                close(connectome_fig);
            end
        end
    end
    
    if plot_tckLengths
        CC_mask = ones(size(vertices.all,1), size(vertices.all,1));
        CC_mask(1:size(vertices.left,1), 1:size(vertices.left,1)) = zeros(size(vertices.left,1), size(vertices.left,1));
        CC_mask(size(vertices.left,1)+1:end, size(vertices.left,1)+1:end) = zeros(size(vertices.right,1), size(vertices.right,1));
        CC_tckLengths_nb = (tckLengths_nb .* CC_mask);
        
        tckLengths_fig = figure('Position', [600 600 1200 600], 'Visible', vis);
        subplot(1,2,1); 
        histogram(tckLengths_nb(tckLengths_nb~=0));
        xlabel('track length (nbr of segments)'); ylabel('counts');
        title('Whole brain');
        
        subplot(1,2,2); 
        histogram(CC_tckLengths_nb(CC_tckLengths_nb~=0));
        xlabel('track length (nbr of segments)'); ylabel('counts');
        title('Inter-hemispheric');
        
        if save_figs
            fname = strcat('tckLengths', save_suffix);
            pause(2)
            export_fig(strcat(PRD, '/connectivity/img_',SUBJ_ID,'/', fname), '-png', '-transparent');
            pause(2);
            export_fig(strcat(PRD, '/connectivity/img_',SUBJ_ID,'/', fname), '-eps');
            pause(2)
            display(strcat(fname, ' figures exported'));
            if close_figs
                close(tckLengths_fig);
            end
        end
    end
    
    if plot_tck_stats    
        tckStats_fig = figure('Position', [600 1200 1200 400]);
        
        subplot(1,3,1);
        histogram(n_isects, -0.5:1:9.5, 'Normalization', 'Probability');
        xlabel('nbr of intersection per bound'); ylabel('Probability');
        
        subplot(1,3,2);
        histogram(stepSz, 0:0.1:5, 'Normalization', 'Probability');
        xlabel('Step size (mm)'); ylabel('Probability');
        
        subplot(1,3,3);
        histogram(tckLengths_mm, 50, 'Normalization', 'Probability');
        xlabel('Track length (mm)'); ylabel('Probability');
        
        if save_figs
            fname = strcat('trackStats_smoothIt', num2str(smoothIter), '_nsamples', num2str(nsamples), '_szBnd', num2str(sz_bound), '_extScaling', num2str(ext_scaling));
            pause(2)
            export_fig(strcat(PRD, '/connectivity/img_',SUBJ_ID,'/', fname), '-png', '-transparent');
            pause(2);
            saveas(tckStats_fig, strcat(PRD, '/connectivity/img_',SUBJ_ID,'/', fname), 'svg'); pause(2);
            display(strcat(fname, ' figure exported'));
            if close_figs
                close(tckStats_fig);
            end
        end
    end
end

% update suffix for saving files
load_suffix = strcat('_', surf_type, '_', connectome_type, '_nsamples', num2str(nsamples), '_smoothIter', num2str(smoothIter), '_szBnd', num2str(sz_bound), '_extScaling', num2str(ext_scaling), '_minTckLen', num2str(minTckLength));
save_suffix = strcat('_',  surf_type, '_', connectome_type, '_nsamples', num2str(nsamples), '_smoothIter', num2str(smoothIter), '_szBnd', num2str(sz_bound), '_extScaling', num2str(ext_scaling), '_minTckLen', num2str(minTckLength),'_10taZsc', num2str(threshold_adjacency*10), trimAdj_all, '_localSpread', num2str(local_spread), '_anisotropy', distAni, num2str(anisotropy*100), '_callosectomy', num2str(callosectomy*100), distCallo, callo_suffix, rand_type, rnd_proba, sep_hemi); %, '_sigmaEigs', num2str(sigma_eigs), '_tolEigs', num2str(tol_eigs));


%%%%%%%%%%%%%%%%%%%%%%%
%% Compute harmonics %%
%%%%%%%%%%%%%%%%%%%%%%%
if compute_harmonics
    if ~exist('connectome', 'var') || reload_connectome
        load(strcat(PRD, '/connectivity/highRes_connectome', load_suffix));
    end
    
    if ~exist('tckLengths_nb', 'var') || reload_connectome
        load(strcat(PRD, '/connectivity/highRes_tckSize', load_suffix));
        if ~exist('tckLengths_nb', 'var')
            tckLengths_nb = tckLengths_nbr;
            clear tckLengths_nbr;
        end
        tckLengths_nb(isnan(tckLengths_nb)) = 0;
        if ~issymmetric(tckLengths_nb)
            tckLengths_nb = (tckLengths_nb + tckLengths_nb')./2;
        end
    end
    
    if ~exist('wm_regmap', 'var')
        % read white matter region mapping
        wm_regmap_lh = load(strcat(PRD, '/surface/lh_white_region_mapping_low.txt'));
        wm_regmap_rh = load(strcat(PRD, '/surface/rh_white_region_mapping_low.txt'));
        wm_regmap = [wm_regmap_lh; wm_regmap_rh];
    end
    
    
    connectome(find(eye(size(connectome)))) = 0;    % remove self connections
    c_idx = find(connectome);                       % indices of connections
    sigma_c = std(connectome(c_idx));               % std dev of connections weights
    mu_c = mean(mean(connectome(c_idx)));           % mean of connections weights
    
    % Adjacency matrix of the connectome
    disp('Trimming connectome...');
	  % weighted trimming (z-score based on weights)
    A_c = zeros(size(connectome));              
    if threshold_adjacency
        zsc = (connectome(c_idx) - mu_c) / sigma_c;
        A_c(c_idx) = zsc > threshold_adjacency;
		% semi-weighted trimming (percentage of nbr of connections (not weight),
		% starting with smallest weights).
    elseif trim_adjacency>0
        A_c = connectome > 0;
        ntrim = round(trim_adjacency * numel(c_idx));
        [vals, vals_idx] = sort(connectome(c_idx), rand_trim_type);
	      connectome(c_idx(vals_idx(1:ntrim))) = 0;
    % track lengths based trimming (percentage of nbr of connections (not weights),
    % starting with longest fibers).
    elseif tckLen_trim_adjacency>0
        ntrim = round(tckLen_trim_adjacency * numel(c_idx));
        [vals, vals_idx] = sort(tckLengths_nb(c_idx), tckLen_trim_type);
        connectome(c_idx(vals_idx(1:ntrim))) = 0;
        A_c = connectome > 0;
    else
        A_c = connectome ~= 0;
    end
    
    if ~issymmetric(A_c)
        A_c = (A_c + A_c') > 0;
    end
    
    % compute local mesh adjacency matrix
    A_local = sparse(size(vertices.all,1),size(vertices.all,1));
    
    % local connectivity : immediate neighbours i.e. local spread = 1
    if local_spread > 0
        A_local_lh = triangulation2adjacency(faces.left); 
        A_local_rh = triangulation2adjacency(faces.right);
        A_local(1:size(vertices.left,1), 1:size(vertices.left,1)) = A_local_lh;
        A_local(size(vertices.left,1)+1:end, size(vertices.left,1)+1:end) = A_local_rh;
        
        % neighbours of neighbours i.e. local spread = 2
        if local_spread > 1
            A_local_2 = sparse(size(vertices.all,1),size(vertices.all,1));
            for i = 1:size(A_local,1)
                ci_idx = find(A_local(i,:));
                for j = ci_idx
                    cj_idx = find(A_local(j,:));
                    A_local_2(i,cj_idx) = 1;
                end
                if ~mod(i, 1000)
                    fprintf('%i%% local spread done. \n\r', round(i*100/size(A_local,1)));
                end
            end
            A_local = (A_local + A_local_2) > 0;
            A_local(find(eye(size(A_local)))) = 0;
        end
    end
    
    % random trimming of local connectivity 
    if anisotropy
      if dist_based_local_trim
        if local_spread==1
          load(fullfile(PRD, 'connectivity', 'R_local_1'));
          R_local = R_local_1;
        elseif local_spread==2
          load(fullfile(PRD, 'connectivity', 'R_local_2'));
          R_local = R_local_2;
        else
          error('choose local_spread to be 1 or 2.')
        end
        c_idx_local = find(A_local>0);
        [rows, cols] = ind2sub(size(A_local), c_idx_local);
        n_lc = numel(c_idx_local);
        n_trim = round(n_lc * anisotropy);      % total number of connections to trim
        local_lengths = R_local(c_idx_local);
        [sorted_local_lengths, local_c_idx] = sort(local_lengths, AniSorting);
        idx_trimmed = local_c_idx(1:n_trim);
        A_local(rows(idx_trimmed), cols(idx_trimmed)) = 0;  % remove connections at trimmed indices
        A_local(cols(idx_trimmed), rows(idx_trimmed)) = 0;  % remove connections at trimmed indices
      else
        c_idx_local = find(A_local>0);
        [rows, cols] = ind2sub(size(A_local), c_idx_local);
        n_lc = numel(c_idx_local);              % total number of local connections
        n_trim = round(n_lc * anisotropy);      % total number of connections to trim
        idx_trimmed = randperm(n_lc, n_trim);   % select randomly indices to be trimmed
        
        A_local(rows(idx_trimmed), cols(idx_trimmed)) = 0;  % remove connections at trimmed indices
        A_local(cols(idx_trimmed), rows(idx_trimmed)) = 0;  % remove connections at trimmed indices
      end
    end
    
    
    % create corpus callosum mask to select inter-henispheric connections
    CC_mask = ones(size(vertices.all,1), size(vertices.all,1));
    CC_mask(1:size(vertices.left,1), 1:size(vertices.left,1)) = zeros(size(vertices.left,1), size(vertices.left,1));
    CC_mask(size(vertices.left,1)+1:end, size(vertices.left,1)+1:end) = zeros(size(vertices.right,1), size(vertices.right,1));

    CC_tckLengths_nb = (tckLengths_nb .* A_c .* CC_mask);

    if plot_CC_stats
        CC_stats_fig = figure('Position', [800 800 1200 400], 'visible', vis);
        subplot(1,3,1);
        histogram(CC_tckLengths_nb(CC_tckLengths_nb~=0));
        xlabel('CC tck lengths (nbr of points)'); ylabel('counts');
        title(sprintf('Before callosectomy n = %i', sum(sum(full(CC_tckLengths_nb~=0)))));
    end

    % trim interhemispheric connections by descending length
    if callosectomy
        CC_idx = find(CC_tckLengths_nb>0);
        CC_pos = CC_tckLengths_nb(CC_idx);
        if dist_based_callo_trim
          [sorted_CC_tckLen, sorted_CC_idx] = sort(CC_pos, CalloSorting);
        else
          sorted_CC_idx = randperm(numel(CC_pos));  % actually not sorted but kept name for compatibility
          sorted_CC_tckLen = CC_pos(sorted_CC_idx); 
        end
        
        trim_sorted_idx = 1:ceil(callosectomy * numel(sorted_CC_tckLen));
        [CC_i, CC_j] = ind2sub(size(A_c), CC_idx(sorted_CC_idx(trim_sorted_idx)));
        A_c(CC_i, CC_j) = 0;
        A_c(CC_j, CC_i) = 0;
    end
    
    % Randomizations
    if randomize_graph_global
        if rand_proba==1
          A_c = randmio_und(A_c, iter_rand);
        else
          n_rows = size(A_c,1); 
          n_rand = floor(rand_proba*n_rows);
          rand_idx = randperm(n_rows, n_rand);
          A_c(rand_idx,rand_idx) = randmio_und(A_c(rand_idx,rand_idx), iter_rand);
        end
    end
    
    if randomize_graph_CC_only
        A_CC_c = A_c(size(vertices.left)+1:end, 1:size(vertices.left,1));
        A_CC_rand = randmio_dir(A_CC_c, iter_rand);
        A_c(size(vertices.left)+1:end, 1:size(vertices.left,1)) = A_CC_rand;
        A_c(1:size(vertices.left,1), size(vertices.left)+1:end) = A_CC_rand';
    end
    
    if randomize_graph_intraHemi_only
        A_left = A_c(1:size(vertices.left), 1:size(vertices.left,1));
        A_left_rand = randmio_und(A_left,iter_rand);
        A_c(1:size(vertices.left), 1:size(vertices.left,1)) = A_left_rand;
        
        A_right = A_c(size(vertices.left)+1:end, size(vertices.left,1)+1:end);
        A_right_rand = randmio_und(A_right,iter_rand);
        A_c(size(vertices.left)+1:end, size(vertices.left,1)+1:end) = A_right_rand;
    end
        
    A = (A_local + A_c) > 0;        % combine local and global adjacency matrices
    
    if plot_CC_stats
        CC_tckLengths_nb = (tckLengths_nb .* A_c .* CC_mask); % update with new A_c (callosal trimmed)
        subplot(1,3,2);
        histogram(CC_tckLengths_nb(CC_tckLengths_nb~=0));
        xlabel('CC tck lengths (nbr of points)'); ylabel('counts');
        title(sprintf('After callosectomy n = %i', sum(sum(full(CC_tckLengths_nb~=0)))));
    end
    
    % degree matrix of adjacency matrix
    D_A = degree(graph(A));
    
    % remove isolated vertices (D_A==0) and subcortical areas ([...] is the array containing the indices of subcortical areas in the region mapping) - /!\ Only Desikan-Killiany atlas supported
    subCtx_arr = ismember(wm_regmap, [0 1 2 3 4 5 6 7 8 9 10 45 46 47 48 49 50 51 52 53]);  
    ctx_idx = find(~subCtx_arr .* ~(D_A==0));
    A_ctx = zeros(size(A));
    A_ctx(ctx_idx, ctx_idx) = A(ctx_idx, ctx_idx);
    
    if plot_CC_stats
        CC_tckLengths_nb = (tckLengths_nb .* A_ctx .* CC_mask); % update with new A_c (callosal trimmed)
        subplot(1,3,3);
        histogram(CC_tckLengths_nb(CC_tckLengths_nb~=0));
        xlabel('CC tck lengths (nbr of points)'); ylabel('counts');
        title(sprintf('After callosectomy \n\r and removing subcortical areas, n = %i', sum(sum(full(CC_tckLengths_nb~=0)))));
        
        % separate figure inter-hemispheric degree projection on cortical surface
        CC_deg_fig = plotNFsurf('CC degree', vertices, faces, degree(graph(CC_tckLengths_nb)), subj_type, vis);
        
        if save_figs
            % save degree
            fname = strcat('CC_tckLengths_nb_degree', save_suffix);
            export_fig(strcat(PRD, '/connectivity/img_',SUBJ_ID,'/', fname), '-png', '-transparent');
            pause(2);
            % save distribution of CC tck lengths
            fname = strcat(PRD, '/connectivity/img_',SUBJ_ID,'/CC_tckLengths_nb_histo', save_suffix);
            pause(1);
            saveas(CC_stats_fig, fname, 'png');
            pause(1);
            saveas(CC_stats_fig, fname, 'svg');
            pause(2);
        end
        if close_figs
            close(CC_stats_fig);
            close(CC_deg_fig);
        end 
    end
    
    % discard isolated clusters disconnected from main graph connectome
    G = graph(A_ctx);
    G_bins = conncomp(G);
    G_idx = find(G_bins==1);        % get node indices of main cluster
   
    if n_graph_chunks == 2
      % find index of second cluster
      ubin = unique(G_bins);
      nbin = zeros(numel(ubin),1);
      for b = 1:numel(ubin)
        nbin(b) = numel(find(G_bins==ubin(b)));
      end
      [nbins, bin_idx] = sort(nbin, 'descend');
      ubin(bin_idx);
      val = ubin(bin_idx(2));
      G_idx2 = find(G_bins==val);
      G_idx = cat(2, G_idx, G_idx2);
    end 

    if plot_graph
        g_fig = figure('Position', [600 600 800 800], 'Visible', vis);
        plot(graph(A_ctx(G_idx, G_idx)));
        if save_figs
            fname = strcat('Graph', save_suffix);
            pause(2);
            export_fig(strcat(PRD, '/connectivity/img_',SUBJ_ID,'/', fname), '-png', '-transparent');
            pause(2);
            g_fig.Renderer = 'Painters';
            saveas(g_fig, strcat(PRD, '/connectivity/img_',SUBJ_ID,'/', fname), 'svg');
            close(g_fig);
        end
    end
    
    if plot_connectomeStats
      plotConnectomeStats(vertices, faces, connectome, A_c, save_figs, threshold_adjacency, PRD, SUBJ_ID, subj_type, save_suffix, vis);
    end

    % extract final degree matrix
    D = diag(degree(graph(A_ctx)));
    
    % compute laplacian
    L = 0.5 * ( (D - A_ctx)' + (D - A_ctx) );

    % perform eigen decomposition
    opts.tol   = tol_eigs; 
    opts.issym = 1;
    opts.disp  = 1;
    [evecs, evals] = eigs(L(G_idx, G_idx), nHarmonics, sigma_eigs, opts);
    
    % reconstruct harmonics in entire surface mesh
    H = zeros(size(connectome,1), nHarmonics);
    H(G_idx, :) = evecs;
  
    % case of harmonics computed for each hemisphere separately
    if separate_hemispheres
      % find index of second cluster
      ubin = unique(G_bins);
      nbin = zeros(numel(ubin),1);
      for b = 1:numel(ubin)
        nbin(b) = numel(find(G_bins==ubin(b)));
      end
      [nbins, bin_idx] = sort(nbin, 'descend');
      ubin(bin_idx);
      val = ubin(bin_idx(2));

      % compute eigenvector-eigenvalue pair of second cluster (i.e. 2nd hemisphere)
      G_idx2 = find(G_bins==val);
      [evecs2, evals2] = eigs(L(G_idx2, G_idx2), nHarmonics, sigma_eigs, opts);
      
      % create harmonics by assembling pairs or eigenvectors using same and different polarity
      H = zeros(size(connectome,1), nHarmonics);
      H(G_idx, 1:2:nHarmonics) = evecs(:,1:(nHarmonics/2));
      H(G_idx, 2:2:nHarmonics) = -evecs(:,1:(nHarmonics/2));
      H(G_idx2,1:2:nHarmonics) = evecs2(:,1:(nHarmonics/2));
      H(G_idx2,2:2:nHarmonics) = evecs2(:,1:(nHarmonics/2));
    end

    % save everything needed for plotting and re-use
    combined = struct('vertices', vertices, 'faces', faces, 'L', sparse(L), 'D_A', D_A, 'eigvals', evals, 'H', H, ...
        'local_vs_global_ratio', sum(sum(A_local(ctx_idx,ctx_idx))) / sum(sum(A_ctx)), 'A_ctx', sparse(A_ctx), 'A_local', sparse(A_local), 'D', sparse(D), 'ta_zsc', threshold_adjacency, ...
        'mu_cc', mu_c, 'sigma_cc', sigma_c, 'ta_w', threshold_adjacency*sigma_c+mu_c);
    if save_combined
        save(strcat(PRD, '/connectivity/combined_harmonics', save_suffix, '.mat'), '-struct', 'combined');
    end
    
    % plot eigenvalues
    fig = figure('Visible', 'off'); bar(diag(evals)); title('Eigenvalues');
    if save_figs
        saveas(fig, strcat(PRD, '/connectivity/img_',SUBJ_ID,'/eigenvalues', save_suffix), 'svg');
    end
    close(fig);
    
    % plot harmonics
    if plot_harmonics
        for h=1:nPlottedHarmonics
            fname = strcat('H', num2str(h, '%02i'), save_suffix);
            hrmnc = plotNFsurf(fname, combined.vertices, combined.faces, combined.H(:,h), subj_type, vis);
            if save_figs
                pause(2);
                export_fig(strcat(PRD, '/connectivity/img_',SUBJ_ID,'/', fname), '-png', '-transparent');
                pause(2);
                close(hrmnc);
                display(strcat(fname, ' image exported'));
            end
        end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Convert High Resolution Harmonics into atlas space %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if convert_into_DSK_atlas
    fprintf('Converting Harmonics from mesh resolution to DSK atlas... \n\r');
    load_suffix = strcat('_', surf_type, '_', connectome_type, '_nsamples', num2str(nsamples), '_smoothIter', num2str(smoothIter), '_szBnd', num2str(sz_bound), '_extScaling', num2str(ext_scaling), '_10taZsc', num2str(threshold_adjacency*10), trimAdj_all, '_localSpread', num2str(local_spread), '_anisotropy', distAni, num2str(anisotropy*100), '_callosectomy', num2str(callosectomy*100), distCallo, callo_suffix,  rand_type, rnd_proba, sep_hemi); %, '_sigmaEigs', num2str(sigma_eigs), '_tolEigs', num2str(tol_eigs));
    save_suffix = strcat('_', surf_type, '_', connectome_type, '_nsamples', num2str(nsamples), '_smoothIter', num2str(smoothIter), '_szBnd', num2str(sz_bound), '_extScaling', num2str(ext_scaling), '_10taZsc', num2str(threshold_adjacency*10), trimAdj_all, '_localSpread', num2str(local_spread), '_anisotropy', distAni, num2str(anisotropy*100), '_callosectomy', num2str(callosectomy*100), distCallo, callo_suffix, rand_type, rnd_proba, sep_hemi); %, '_10thZscMI', num2str(thZscMI*10)); %, '_sigmaEigs', num2str(sigma_eigs), '_tolEigs', num2str(tol_eigs));
    if (~exist('combined', 'var') || reload_harmonics)
        combined = load(strcat(PRD, '/connectivity/combined_harmonics', load_suffix, '.mat'));
    end
    
    if ~exist('wm_regmap', 'var')
        % read white matter region mapping
        wm_regmap_lh = load(strcat(PRD, '/surface/lh_white_region_mapping_low.txt'));
        wm_regmap_rh = load(strcat(PRD, '/surface/rh_white_region_mapping_low.txt'));
        wm_regmap = [wm_regmap_lh; wm_regmap_rh];
    end
    
    ROIs = unique(wm_regmap);
    nROI = numel(ROIs);
    H_DSK = zeros(nHarmonics,nROI);
    for h = 1:nHarmonics
        for i = 1:nROI
            H_DSK(h,i) = mean(combined.H(find(wm_regmap==ROIs(i)), h));
        end
    end
    
    if plot_converted_harmonics
        for h=1:nPlottedHarmonics
            tmpH = zeros(1,size(combined.H,1));
            for i = 1:nROI
                tmpH(find(wm_regmap==ROIs(i))) = H_DSK(h,i);
            end
            fname = strcat('H_DSK', num2str(h, '%02i'), save_suffix);
            hrmnc = plotNFsurf(fname, combined.vertices, combined.faces, tmpH', subj_type, vis);
            if save_figs
                pause(2);
                export_fig(strcat(PRD, '/connectivity/img_',SUBJ_ID,'/', fname), '-png', '-transparent');
                pause(2);
                hrmnc.Renderer = 'Painters';
                saveas(hrmnc, strcat(PRD, '/connectivity/img_', SUBJ_ID,'/', fname), 'svg'); pause(1)
                close(hrmnc);
                display(strcat(fname, ' image exported'));
            end
        end
    end
    
    if re_save_combined
        combined.H_DSK = sparse(H_DSK);
        if isfield(combined, 'A_ctx')
            if ~issparse(combined.A_ctx)
                combined.A_ctx = sparse(combined.A_ctx);
            end
        end
        if isfield(combined, 'A')
            if ~issparse(combined.A)
                combined.A = sparse(combined.A);
            end
        end
        if isfield(combined,'D')
            if ~issparse(combined.D)
                combined.D = sparse(combined.D);
            end
        end
        if isfield(combined,'L')
            if ~issparse(combined.L)
                combined.L = sparse(combined.L);
            end
        end
        save(strcat(PRD, '/connectivity/combined_harmonics', save_suffix, '.mat'), '-struct', 'combined');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mutual information with RSNs based on region mapping %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if MI_RSN
    fprintf('Computing Mutual Information... \n\r');
    load_suffix = strcat('_', surf_type, '_', connectome_type, '_nsamples', num2str(nsamples), '_smoothIter', num2str(smoothIter), '_szBnd', num2str(sz_bound), '_extScaling', num2str(ext_scaling), '_minTckLen', num2str(minTckLength), '_10taZsc', num2str(threshold_adjacency*10), trimAdj_all, '_localSpread', num2str(local_spread), '_anisotropy', distAni, num2str(anisotropy*100), '_callosectomy', num2str(callosectomy*100), distCallo, callo_suffix, rand_type, rnd_proba, sep_hemi); %, '_sigmaEigs', num2str(sigma_eigs), '_tolEigs', num2str(tol_eigs));
    save_suffix = strcat('_', surf_type, '_', connectome_type, '_nsamples', num2str(nsamples), '_smoothIter', num2str(smoothIter), '_szBnd', num2str(sz_bound), '_extScaling', num2str(ext_scaling), '_minTckLen', num2str(minTckLength), '_10taZsc', num2str(threshold_adjacency*10), trimAdj_all, '_localSpread', num2str(local_spread), '_anisotropy', distAni, num2str(anisotropy*100), '_callosectomy', num2str(callosectomy*100), distCallo, callo_suffix, rand_type, rnd_proba, sep_hemi, '_10thZscMI', num2str(thZscMI*10)); %, '_sigmaEigs', num2str(sigma_eigs), '_tolEigs', num2str(tol_eigs));
    if (~exist('combined', 'var') || reload_harmonics)
        combined = load(strcat(PRD, '/connectivity/combined_harmonics', load_suffix, '.mat'));
    end
    
    if ~exist('wm_regmap', 'var')
        % read white matter region mapping
        wm_regmap_lh = load(strcat(PRD, '/surface/lh_white_region_mapping_low.txt'));
        wm_regmap_rh = load(strcat(PRD, '/surface/rh_white_region_mapping_low.txt'));
        wm_regmap = [wm_regmap_lh; wm_regmap_rh];
    end
    
    % select DMN areas ([...] is the array containing the indices of DMN areas in the region mapping: isthmus cing, mPFC, middle temp, PCC, angular gyrus)
    RSN_version = '2020_';  %'2018_' <-- older versions of the code used less regions in DMN
     
    % region wise
    if strcmp(RSN, 'DMN')
      if strcmp(RSN_version, '2020_')
        rsn_lh = [13 17 19 24 32 34 35 37 39]; 
        rsn_rh = [60 62 67 75 77 78 80 82];
      else
        rsn_lh = [19 23 24 32 40]; 
        rsn_rh = [62 66 67 75 83];
      end
    elseif strcmp(RSN, 'SM')
      rsn_lh = [31 33 26];
      rsn_rh = [74 76 69];
    elseif strcmp(RSN, 'FP')
      rsn_lh = [12 16 18 36];
      rsn_rh = [55 59 61 79];
    elseif strcmp(RSN, 'Vis')
      rsn_lh = [20 30 14 22];
      rsn_rh = [63 73 57 65];
    elseif strcmp(RSN, 'L')
      rsn_lh = [15 21 23 41 42];
      rsn_rh = [58 64 66 84 85];
    elseif strcmp(RSN, 'VA')
      rsn_lh = [27 44];
      rsn_rh = [70 87];
    elseif strcmp(RSN, 'DA')
      rsn_lh = [38 13];
      rsn_rh = [81 56];
    end
    
    % based on Yeo et al. 2011
    if RSN_Yeo2011
      if ~contains(SUBJ_ID, 'fsaverage')
        error("Yeo2011 parcellaltion only available for fsaverage surfaces");
      end
      [v.left, l.left, c.left] = read_annotation(strcat(PRD,'/label/lh.Yeo2011_7Networks_N1000.annot'));
      [v.right, l.right, c.right] = read_annotation(strcat(PRD,'/label/rh.Yeo2011_7Networks_N1000.annot'));
      wm_regmap_lh = l.left;
      wm_regmap_rh = l.right;
      wm_regmap = [wm_regmap_lh; wm_regmap_rh];
      if strcmp(RSN, 'DA') 
        rsn_lh = [947712];
        rsn_rh = [941712];
      elseif strcmp(RSN, 'FP') 
        rsn_lh = [2266342];
        rsn_rh = [2266342];
      elseif strcmp(RSN, 'DMN')
        rsn_lh = [5127885];
        rsn_rh = [5127885];
      elseif strcmp(RSN, 'Vis')
        rsn_lh = [8786552];
        rsn_rh = [8786552];
      elseif strcmp(RSN, 'L')
        rsn_lh = [10811612];
        rsn_rh = [10811612];
      elseif strcmp(RSN, 'SM') 
        rsn_lh = [11829830];
        rsn_rh = [11829830];
      elseif strcmp(RSN, 'VA')
        rsn_lh = [16399044];
        rsn_rh = [16399044];
      end 
      RSN_version = '_Yeo2011';
    end

    rsn_lh_arr = ismember(wm_regmap_lh, rsn_lh);
    rsn_rh_arr = ismember(wm_regmap_rh, rsn_rh);
    rsn_lh_idx = find(rsn_lh_arr);
    rsn_rh_idx = find(rsn_rh_arr);
   
    % whole brain
    wm_rsn_arr = ismember(wm_regmap, cat(2, rsn_lh, rsn_rh));  
    wm_rsn_idx = find(wm_rsn_arr);
    
    
    rsn_arr = [rsn_lh_arr; rsn_rh_arr; wm_rsn_arr];
    lh_idxs = 1:numel(wm_regmap_lh);
    rh_idxs = numel(wm_regmap_lh)+1:numel(wm_regmap);
    all_idxs = 1:numel(wm_regmap);
    idxs = {lh_idxs, rh_idxs, all_idxs};
    
    if plot_RSN_mask
        RSN_fig = plotNFsurf(RSN, combined.vertices, combined.faces, double(wm_rsn_arr), 'HCP', vis);
        if save_figs
            pause(1)
            print(RSN_fig, strcat(PRD, '/connectivity/img_', SUBJ_ID,'/', rsn_type, save_suffix), '-dpng'); pause(1)
        end
        if strcmp(vis, 'off')
            close(RSN_fig);
        end
    end

    
    % Compute MI with DMN using different methods
    for m = 1:numel(nBinsMI)
        % transform harmonics into binary vectors (using z-score thresholding at value thZscMI)
        if nBinsMI(m)==1
            for h=1:nHarmonics
                for hemi = 1:3  % left, right, both
                    zharm = ( combined.H(:,h) - mean(combined.H(:,h)) ) / std(combined.H(:,h));
                    zharm_bin  = zharm > thZscMI;

                    J = [ sum(~wm_rsn_arr(idxs{hemi}) & ~zharm_bin(idxs{hemi})), sum(~wm_rsn_arr(idxs{hemi}) & zharm_bin(idxs{hemi})); ...
                          sum( wm_rsn_arr(idxs{hemi}) & ~zharm_bin(idxs{hemi})), sum( wm_rsn_arr(idxs{hemi}) & zharm_bin(idxs{hemi})) ] / numel(zharm_bin(idxs{hemi}));
                    JJ = sum( sum( J .* log2(J./(sum(J,2)*sum(J,1))) ) );

                    % sign reversal 
                    zharm_bin  = -zharm > thZscMI;
                    K = [ sum(~wm_rsn_arr(idxs{hemi}) & ~zharm_bin(idxs{hemi})), sum(~wm_rsn_arr(idxs{hemi}) & zharm_bin(idxs{hemi})); ...
                          sum( wm_rsn_arr(idxs{hemi}) & ~zharm_bin(idxs{hemi})), sum( wm_rsn_arr(idxs{hemi}) & zharm_bin(idxs{hemi})) ] / numel(zharm_bin(idxs{hemi}));
                    KK = sum( sum( K .* log2(K./(sum(K,2)*sum(K,1))) ) );
                    
                    MI(m,hemi,h) = max(JJ,KK);

                    CC(h) = corr(double(wm_rsn_arr), double(zharm_bin));
                end
            end
            MI_suffix = 'v3';
        else
            % discretize harmonics into n bins
            for h = 1:nHarmonics
                for hemi = 1:3  % left, right, both
                    [nh,eh] = histcounts(combined.H(:,h), nBinsMI(m), 'Normalization', 'probability');
                    [nd,ed] = histcounts(wm_rsn_arr, 2, 'Normalization', 'probability');
                    for i=1:nBinsMI(m)  % binarized (i.e. more "continuous") data (Harmonics)
                        for j=1:2       % binary data (DMN)
                            jdp(i,j) = sum((combined.H(idxs{hemi},h)>eh(i)) .* (combined.H(idxs{hemi},h)<=eh(i+1)) .* (wm_rsn_arr(idxs{hemi})>ed(j)) .* (wm_rsn_arr(idxs{hemi})<=ed(j+1)))/numel(wm_rsn_arr(idxs{hemi}));
                            idp(i,j) = nh(i)*nd(j);
                        end
                    end
                    idx = jdp~=0;
                    MI(m,hemi,h) = sum( sum( jdp(idx) .* log2(jdp(idx) ./ idp(idx)) ) );
                end
            end
            MI_suffix = strcat('nBins', num2str(nBinsMI(m)));
        end
        MI_m = permute(MI(m,:,:), [3 2 1]);
        MI_fig = figure('Position', [0 1000 3200 300], 'Visible', vis); 
        bar(1:nHarmonics, MI_m);
        ylim(MI_ylim);
        xlabel('Harmonics'); ylabel('Mutual Information');
        legend({'left hemi', 'right hemi', 'both'})
        if save_figs
            pause(1)
            fname = strcat(PRD, '/connectivity/img_', SUBJ_ID,'/MI_', rsn_type, MI_suffix, save_suffix);
            print(MI_fig, fname, '-depsc');
            print(MI_fig, fname, '-dpng');
            print(MI_fig, fname, '-dsvg');
            
        end
        if save_MI
            save(strcat(PRD, '/connectivity/MI_',rsn_type, RSN_version, MI_suffix, save_suffix), 'MI_m');
        end
        if strcmp(vis, 'off')
            close(MI_fig);
        end
    end
end


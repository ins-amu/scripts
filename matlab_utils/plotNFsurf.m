% Creates a visualization of a field (2D, 3D) on the cortical surface described by vertices and faces (see
% freesurfer)
%
% dependencies: toolbox_graph, subtightplot
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

function h = plotNFsurf(fname, vertices, faces, field, subj_type, visible)

subplot = @(m,n,p) subtightplot (m, n, p, [0.02 0.05], [0.02 0.02], [0.05 0.05]);

%inds.left = 1:length(vertices.left);
%inds.right = length(vertices.left)+1:length(vertices.all);  

%% preprocessing
% make sure data is zero mean (for harmonics)
%field = field - mean(field);

% get colormap
n_colors = 105; % discretization of colormap
if      strcmp(subj_type,'HD')
            cbrewer_cmap = cbrewer('div','BrBG',n_colors);

elseif  strcmp(subj_type,'HCP')
            cbrewer_cmap = cbrewer('div','RdBu',n_colors);

elseif  strcmp(subj_type,'RegMap')
            cbrewer_cmap = cbrewer('qual','Paired', 86);
    
else
            cbrewer_cmap = cbrewer('div','PRGn',n_colors);
end
    

if size(field,2)==1
    options.face_vertex_color = field;
else
    options.normal = field;
end

h = figure( 'Name', fname,   'Units', 'Normalized', ...
            'Position', [0.3 0.3 0.8 0.4], 'Visible', visible );
cmax = max(field);
cmin = min(field);
c = max(abs([cmax, cmin]));

% Left-Right Hemisphere boundary index, assumed exact same size of hemispheres for now
LRH_boundary = round(numel(field)/2);  

% left hemisphere - lateral view
subplot(1,5,1);
options.face_vertex_color = field(1:LRH_boundary);
plot_mesh(vertices.left, faces.left, options);
caxis([-c c]); 
%caxis([cmin cmax]); 
colormap(cbrewer_cmap); view([-90 0]); camlight;

% left hemisphere - medial view
subplot(1,5,2);
options.face_vertex_color = field(1:LRH_boundary);
plot_mesh(vertices.left, faces.left, options);
caxis([-c c]); 
%caxis([cmin cmax]);
colormap(cbrewer_cmap); view([90 0]); camlight;

% top view
subplot(1,5,3);
options.face_vertex_color = field(1:LRH_boundary);
plot_mesh(vertices.left, faces.left, options);
caxis([-c c]); 
%caxis([cmin cmax]);
hold on;
options.face_vertex_color = field(LRH_boundary+1:end);
plot_mesh(vertices.right, faces.right, options);
caxis([-c c]); 
%caxis([cmin cmax]);
colormap(cbrewer_cmap); camlight;

% right hemisphere - medial view
subplot(1,5,4);
options.face_vertex_color = field(LRH_boundary+1:end);
plot_mesh(vertices.right, faces.right, options);
caxis([-c c]); 
%caxis([cmin cmax]);
colormap(cbrewer_cmap); view([-90 0]); camlight;

% right hemisphere - lateral view
subplot(1,5,5);
options.face_vertex_color = field(LRH_boundary+1:end);
plot_mesh(vertices.right, faces.right, options);
caxis([-c c]); 
%caxis([cmin cmax]);
colormap(cbrewer_cmap); view([90 0]); camlight;

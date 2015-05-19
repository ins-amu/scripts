FREESURFER_HOME = fullfile(getenv('FREESURFER_HOME'), '/')
FS = fullfile(getenv('FS'), '/')
PRD = fullfile(getenv('PRD'), '/')

if (~isdeployed)
    addpath(fullfile(FREESURFER_HOME, 'matlab'))
end

if ~exist('rl', 'var')
    if ~exist(fullfile(PRD, 'surface', 'rh.pial.asc'))
        rl='lh'
    else
        rl='rh'
    end
end

corr_right = load([rl, '_ref_table.txt']);
SUBJ_ID = getenv('SUBJ_ID')
[v, L, ct] = read_annotation(fullfile(FS,SUBJ_ID, 'label', [rl, '.aparc.annot']));
cd(PRD)
a = load(fullfile('surface', [rl ,'_vertices_low.txt']));
b = load(fullfile('surface', [rl, '_triangles_low.txt']));
c = load(fullfile('surface', [rl, '_vertices_high.txt']));
reg_map = zeros(size(a,1),1);
not_found = [];
for i=1:size(a,1)
    i;
    [g,e] = min(abs(c(:,1)-a(i,1))+abs(c(:,2)-a(i,2))+abs(c(:,3)-a(i,3)));
    find_tab = find(corr_right(:,6) == L(e));
    if isempty(find_tab)
        not_found = [i, not_found];
    else
        reg_map(i) = corr_right(find_tab,5);
    end
end
not_found
save(fullfile('surface', [rl, '_region_mapping_low_not_corrected.txt']),'reg_map', '-ascii' );

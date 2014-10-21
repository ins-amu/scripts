FS = getenv('FS')

if (~isdeployed)
    addpath([FS, '/matlab'])
end

if ~exist('rl', 'var')
    if ~exist([PRD, '/surface/', 'rh.pial.asc'])
        rl='lh'
    else
        rl='rh'
    end
end

PRD = getenv('PRD')
corr_right = load([rl, '_ref_table.txt']);
SUBJ_ID = getenv('SUBJ_ID')
[v, L, ct] = read_annotation([FS,'/',SUBJ_ID, '/label/', rl, '.aparc.annot']);
cd(PRD)
a = load(['surface/',rl ,'_vertices_low.txt']);
b = load(['surface/', rl, '_triangles_low.txt']);
c = load(['surface/', rl, '_vertices_high.txt']);
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
save(['surface/', rl, '_region_mapping_low_not_corrected.txt'],'reg_map', '-ascii' );

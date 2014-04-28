PRD = getenv('PRD')
cd(PRD)
FS = getenv('FS')
SUBJ_ID = getenv('SUBJ_ID')
[v, L, ct] = read_annotation([FS,'/',SUBJ_ID, '/label/lh.aparc.annot']); 
a = load('surface/lh_vertices_low.txt');
b = load('surface/lh_triangles_low.txt');
c = load('surface/lh_vertices_high.txt');
corr_left = load('scripts/left_hemi_ref_table.txt');
reg_map = zeros(size(a,1),1);
not_found = [];
for i=1:size(a,1)
	i;
	[g,e] = min(abs(c(:,1)-a(i,1))+abs(c(:,2)-a(i,2))+abs(c(:,3)-a(i,3)));
	find_tab = find(corr_left(:,6) == L(e));
	if isempty(find_tab)
		not_found = [i, not_found];
	else
		reg_map(i) = corr_left(find_tab,5);
	end
end
not_found

save('surface/lh_region_mapping_low_not_corrected.txt','reg_map', '-ascii' );

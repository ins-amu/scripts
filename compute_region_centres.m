if (~isdeployed)
addpath('read_and_write_func');
end
PRD = getenv('PRD')
SUBJ_ID = getenv('SUBJ_ID')
parcel = getenv('parcel')

g = load_nii([PRD, '/connectivity_regions/region_parcellation.nii']);
corr_mat = load(['parcellations/correspondance_mat_', parcel, '.txt']);
data = g.img;

list_region = unique(data);
centres = zeros(size(list_region, 1)-1, 4);

for j=2:size(list_region, 1)
    list_region(j)
    ind = corr_mat(find(corr_mat(:, 1)==list_region(j)), 2)
    [a, b, c] = ind2sub(size(data), find(data==list_region(j)));
    centres(ind, 2:4) = [mean(a), mean(b), mean(c)];
    centres(ind, 1) = int32(list_region(j));
end

% save([PRD, '/', SUBJ_ID, '_regions/connectivity/centres.txt'], 'centres', '-ascii');
fid = fopen([PRD, '/', SUBJ_ID, '_regions/connectivity/centres.txt'], 'w');
fprintf(fid, '%d %.3f %.3f %.3f\n', centres');
fclose(fid);  

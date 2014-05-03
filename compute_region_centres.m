if (~isdeployed)
addpath('read_and_write_func');
end
PRD = getenv('PRD')
SUBJ_ID = getenv('SUBJ_ID')
parcel = getenv('parcel')

g = load_nii([PRD, '/connectivity_regions/region_parcellation.nii'], [], [], [], [], [], 0.5);
% h = load_untouch_nii([PRD, '/connectivity_regions/region_parcellation.nii']);
corr_mat = load(['parcellations/correspondance_mat_', parcel, '.txt']);
data = g.img;
% datah = h.img;

% affine = [h.hdr.hist.srow_x; h.hdr.hist.srow_y; h.hdr.hist.srow_z]
% rot = affine(1:3, 1:3); 
% dime = h.hdr.dime.dim(2:4);
% trans = zeros(3,1);
% for i = 1:3
%     [k,j] = max(abs(rot(i,:)));
%     rot(abs(rot)<0.1) = 0;
%     if rot(i, j)<0
%         trans(i) = dime(j);
%     end
% end

list_region = corr_mat(:, 1);
centres = zeros(size(list_region, 1)-1, 4);
% centresh = zeros(size(list_region, 1)-1, 4);

for j=1:size(list_region, 1)
    list_region(j)
    [a, b, c] = ind2sub(size(data), find(data==list_region(j)));
    centres(j, 2:4) = [mean(a), mean(b), mean(c)];
    centres(j, 1) = int32(list_region(j));
end

% for j=1:size(list_region, 1)
%     list_region(j)
%     [a, b, c] = ind2sub(size(datah), find(datah==list_region(j)));
%     centresh(j, 2:4) = rot*[mean(a), mean(b), mean(c)]'+trans;
%     centresh(j, 1) = int32(list_region(j));
% end

% save([PRD, '/', SUBJ_ID, '_regions/connectivity/centres.txt'], 'centres', '-ascii');
fid = fopen([PRD, '/', SUBJ_ID, '_regions/connectivity/centres.txt'], 'w');
fprintf(fid, '%d %.3f %.3f %.3f\n', centres');
fclose(fid);  

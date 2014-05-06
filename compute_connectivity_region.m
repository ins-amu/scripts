if (~isdeployed)
addpath('read_and_write_func')
end
PRD = getenv('PRD')
SUBJ_ID = getenv('SUBJ_ID')
number_tracks = str2num(getenv('number_tracks'))
parcel = getenv('parcel')

g = load_untouch_nii([PRD, '/connectivity_regions/region_parcellation_2_diff.nii.gz']);
data = g.img;
size_img = size(g.img);
corr_mat = load(['parcellations/correspondance_mat_', parcel , '.txt']);
list_region = unique(data);
size_parcel = size(corr_mat, 1) 
res = zeros(size_parcel, size_parcel);
res_length = zeros(size_parcel,size_parcel);
num_tracks = 100000;
affin = [g.hdr.hist.srow_x; g.hdr.hist.srow_y; g.hdr.hist.srow_z]
r = inv(affin(1:3,1:3));
j=0;
for nt =1:number_tracks
'iteration'
nt
tracks = read_mrtrix_tracks(sprintf([PRD, '/connectivity_regions/whole_brain_%d.tck'],nt));
for i=1:num_tracks
a = tracks.data{i};
c = -affin(1:3,4);
d= repmat(c',size(a,1),1);
e = abs(r*(a+d)')+1;
%start_point = data(round(e(1,1)),round(size_img(2)-e(2,1)),round(e(3,1)));
%end_point =  data(round(e(1,end)),round(size_img(2)-e(2,end)),round(e(3,end)));
%start_point = data(round(e(1,1)),round(e(2,1)),round(e(3,1)));
%end_point =  data(round(e(1,end)),round(e(2,end)),round(e(3,end)));
countei = 0;
start_point=0;
while start_point==0 && countei < 10
countei = countei + 1;
start_point = data(round(e(1,countei)),round(e(2,countei)),round(e(3,countei)));
end
countei = size(e, 2);
end_point=0;
while end_point==0 && countei > size(e, 2) -10
countei = countei - 1;
end_point = data(round(e(1,countei)),round(e(2,countei)),round(e(3,countei)));
end
if start_point >0 & end_point > 0
start_ind = corr_mat(find(corr_mat(:,1)==start_point),2);
end_ind = corr_mat(find(corr_mat(:,1)==end_point),2);
if start_ind >0 & end_ind > 0
j = j+1;
res(start_ind, end_ind) = 1 + res(start_ind, end_ind);
res_length(start_ind, end_ind) = size(e, 2) + res_length(start_ind, end_ind);
end
end
end
end

'number of tracts'
j
save([PRD, '/connectivity_regions/raw_connectivity.mat'], 'res');
% postprocessing
%distance between voxels in mm: 0.04mm
res_length = 0.04 .* res_length;
% to have the right count
length_mat = res_length./res;
length_mat(isnan(length_mat))=0;
% to avoid the biais toward longer tracts
connectivity_mat =  res./length_mat;
connectivity_mat(isnan(connectivity_mat)) = 0;
% eliminate diagonal
for i=1:size_parcel
connectivity_mat(i,i) = 0;
end

%figure(1)
%imshow(res./max(max(res)), 'Colormap', jet(25))
%figure()
%imshow(log(res)./max(max(log(res))), 'Colormap', jet(25))
%figure(3)
%imshow((res_length./res)./max(max((res_length./res))), 'Colormap', jet(25))
%figure()
%imshow(log(res_length./res)./max(max(log(res_length./res))), 'Colormap', jet(25))
% to get the average length
%length_mat = res_length./res;
% to compensate for the biais in favor of longer fibers
%connectivity_mat =  res./length_mat;
%connectivity_mat = res
%connectivity_mat(isnan(connectivity_mat)) = 0;
%length_mat(isnan(length_mat))=0;
% f1 = figure()
% imshow(log(length_mat)./max(max(log(length_mat))), 'Colormap', jet(255))
% f2 = figure()
% imshow(log(connectivity_mat)./max(max(log(connectivity_mat))), 'Colormap', jet(255))
% saveas(f1,[PRD, '/connectivity/length_1.jpg'],'jpg')
% saveas(f2,[PRD, '/connectivity/connectivity_1.jpg'],'jpg')

save([PRD, '/', SUBJ_ID, '_regions/connectivity/weights.txt'], 'connectivity_mat', '-ascii')
save([PRD, '/', SUBJ_ID, '_regions/connectivity/tract_lengths.txt'], 'length_mat', '-ascii')

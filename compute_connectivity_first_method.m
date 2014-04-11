addpath('read_and_write_func')
PRD = getenv('PRD')
SUBJ_ID = getenv('SUBJ_ID')
res = zeros(88,88);
res_length = zeros(88,88);

g = load_untouch_nii([PRD, '/connectivity/aparcaseg_2_diff.nii.gz']);
size_img = size(g.img);
num_tracks = 100000;
affin = [g.hdr.hist.srow_x; g.hdr.hist.srow_y; g.hdr.hist.srow_z]
r = inv(affin(1:3,1:3));
data = g.img;
corr_mat = load('correspondance_mat.txt');
j=0;
for nt =1:3
'iteration'
nt
tracks = read_mrtrix_tracks(sprintf([PRD, '/connectivity/whole_brain_%d.tck'],nt));
for i=1:num_tracks
a = tracks.data{i};
c = -affin(1:3,4);
d= repmat(c',size(a,1),1);
e = abs(r*(a+d)')+1;
start_point = data(round(e(1,1)),round(size_img(2)-e(2,1)),round(e(3,1)));
end_point =  data(round(e(1,end)),round(size_img(2)-e(2,end)),round(e(3,end)));
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

%figure(1)
%imshow(res./max(max(res)), 'Colormap', jet(25))
%figure()
%imshow(log(res)./max(max(log(res))), 'Colormap', jet(25))
%figure(3)
%imshow((res_length./res)./max(max((res_length./res))), 'Colormap', jet(25))
%figure()
%imshow(log(res_length./res)./max(max(log(res_length./res))), 'Colormap', jet(25))
% to get the average length
length_mat = res_length./res;
% to compensate for the biais in favor of longer fibers
connectivity_mat =  res./length_mat;
connectivity_mat(isnan(connectivity_mat)) = 0;
length_mat(isnan(length_mat))=0;
% f1 = figure()
% imshow(log(length_mat)./max(max(log(length_mat))), 'Colormap', jet(255))
% f2 = figure()
% imshow(log(connectivity_mat)./max(max(log(connectivity_mat))), 'Colormap', jet(255))
% saveas(f1,[PRD, '/connectivity/length_1.jpg'],'jpg')
% saveas(f2,[PRD, '/connectivity/connectivity_1.jpg'],'jpg')
save([PRD, '/', SUBJ_ID, '/connectivity/weights_method1.txt'], 'connectivity_mat', '-ascii')
save([PRD, '/', SUBJ_ID, '/connectivity/tracts_method1.txt'], 'length_mat', '-ascii')

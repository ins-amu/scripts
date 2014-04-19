addpath('read_and_write_func/')
PRD = getenv('PRD')
SUBJ_ID = getenv('SUBJ_ID')
res = zeros(88,88);
res_length = zeros(88,88);
not_found = 0;
num_tracks = 100000;
g = load_untouch_nii([PRD, '/connectivity/aparcaseg_2_diff.nii.gz']);
size_img = size(g.img);
affin = [g.hdr.hist.srow_x; g.hdr.hist.srow_y; g.hdr.hist.srow_z]
r = inv(affin(1:3,1:3));
data = g.img;
corr_mat = load('correspondance_mat.txt');
j=0
for nt =1:10
'iteration'
nt
tracks = read_mrtrix_tracks(sprintf([PRD, '/connectivity/whole_brain_%d.tck'],nt));
    for i=1:num_tracks
    ind = [];
    uind = [];
    point = [];
    a = tracks.data{i};
    c = -affin(1:3,4);
    d= repmat(c',size(a,1),1);
    e = abs(r*(a+d)')+1;
        for i = 1:size(e,2)
        %point = [point, data(round(e(1,i)),round(size_img(2)-e(2,i)),round(e(3,i)))];
        point = [point, data(round(e(1,i)),round(e(2,i)),round(e(3,i)))];
		end
        [upoint, indpoint] = unique(point);
            if size(upoint,2) > 2
                for k=2:size(upoint,2)
                curr = corr_mat(find(corr_mat(:,1)==upoint(k)),2);
                    if ~isempty(curr)
                    ind = [ind, curr];
                    uind = [uind, indpoint(k)];
                    end
                end
            end


        if size(ind,2)>1
            j = j+1;
            dist_curr = abs(repmat(uind',1, size(uind,2)) - repmat(uind, size(uind,2),1));
            res_length(ind, ind)=dist_curr + res(ind, ind);
            res(ind, ind) = 1 + res(ind, ind);
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
%to get the average length
length_mat = res_length./res;
%to compensate for the biais in favor of longer fibers
%connectivity_mat =  res./length_mat;
connectivity_mat = res
connectivity_mat(isnan(connectivity_mat)) = 0;
length_mat(isnan(length_mat))=0;
% f1 = figure()
% imshow(log(length_mat)./max(max(log(length_mat))), 'Colormap', jet(255))
% f2 = figure()
% imshow(log(connectivity_mat)./max(max(log(connectivity_mat))), 'Colormap', jet(255))
% saveas(f1,[PRD, '/connectivity/length_2.jpg'],'jpg')
% saveas(f2,[PRD, '/connectivity/connectivity_2.jpg'],'jpg')
save([PRD, '/', SUBJ_ID, '/connectivity/weights_method3.txt'], 'connectivity_mat', '-ascii')
save([PRD, '/', SUBJ_ID, '/connectivity/tracts_method3.txt'], 'length_mat', '-ascii')

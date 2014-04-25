addpath('read_and_write_func');
PRD = getenv('PRD')
SUBJ_ID = getenv('SUBJ_ID')

g = load_nii([PRD, '/connectivity/region_parcellation.nii']);
data = g.img;
list_region = unique(data);
mean_centres = zeros(size(list_region, 1)-1,3);
for j=2:size(list_region, 1)
    list_region(j)
    [a, b, c] = ind2sub(size(data), find(data==list_region(j)));
    mean_centres(j-1, :) = [mean(a), mean(b), mean(c)];
end

centres = [double(list_region), mean_centres]
save([PRD, SUBJ_ID, '_regions/centres.txt'], 'centres', '-ascii') 

if (~isdeployed) 

addpath('read_and_write_func')   
addpath('uniform_parcellate')

end

curr_K = getenv('curr_K')
curr_K = str2num(curr_K);
K = log(curr_K)/log(2) 
PRD = fullfile(getenv('PRD'), '/')
SUBJ_ID = getenv('SUBJ_ID')

% get rid of subcortical regions and white matter

mask_unprocessed = load_untouch_nii(fullfile(PRD, 'connectivity', 'aparcaseg_2_diff.nii.gz')); 
dat = mask_unprocessed.img;
dat(ind2sub(size(dat), find(dat ==   0)))=0;
dat(ind2sub(size(dat), find(dat ==   2)))=0;
dat(ind2sub(size(dat), find(dat ==   4)))=0;
dat(ind2sub(size(dat), find(dat ==   5)))=0;
dat(ind2sub(size(dat), find(dat ==   7)))=0;
dat(ind2sub(size(dat), find(dat ==   8)))=0;
dat(ind2sub(size(dat), find(dat ==  10)))=0;
dat(ind2sub(size(dat), find(dat ==  11)))=0;
dat(ind2sub(size(dat), find(dat ==  12)))=0;
dat(ind2sub(size(dat), find(dat ==  13)))=0;
dat(ind2sub(size(dat), find(dat ==  14)))=0;
dat(ind2sub(size(dat), find(dat ==  15)))=0;
dat(ind2sub(size(dat), find(dat ==  16)))=0;
dat(ind2sub(size(dat), find(dat ==  17)))=0;
dat(ind2sub(size(dat), find(dat ==  18)))=0;
dat(ind2sub(size(dat), find(dat ==  24)))=0;
dat(ind2sub(size(dat), find(dat ==  26)))=0;
dat(ind2sub(size(dat), find(dat ==  28)))=0;
dat(ind2sub(size(dat), find(dat ==  30)))=0;
dat(ind2sub(size(dat), find(dat ==  31)))=0;
dat(ind2sub(size(dat), find(dat ==  41)))=0;
dat(ind2sub(size(dat), find(dat ==  43)))=0;
dat(ind2sub(size(dat), find(dat ==  44)))=0;
dat(ind2sub(size(dat), find(dat ==  46)))=0;
dat(ind2sub(size(dat), find(dat ==  47)))=0;
dat(ind2sub(size(dat), find(dat ==  49)))=0;
dat(ind2sub(size(dat), find(dat ==  50)))=0;
dat(ind2sub(size(dat), find(dat ==  51)))=0;
dat(ind2sub(size(dat), find(dat ==  52)))=0;
dat(ind2sub(size(dat), find(dat ==  53)))=0;
dat(ind2sub(size(dat), find(dat ==  54)))=0;
dat(ind2sub(size(dat), find(dat ==  58)))=0;
dat(ind2sub(size(dat), find(dat ==  60)))=0;
dat(ind2sub(size(dat), find(dat ==  62)))=0;
dat(ind2sub(size(dat), find(dat ==  63)))=0;
dat(ind2sub(size(dat), find(dat ==  72)))=0;
dat(ind2sub(size(dat), find(dat ==  77)))=0;
dat(ind2sub(size(dat), find(dat ==  80)))=0;
dat(ind2sub(size(dat), find(dat ==  85)))=0;
dat(ind2sub(size(dat), find(dat == 251)))=0;
dat(ind2sub(size(dat), find(dat == 252)))=0;
dat(ind2sub(size(dat), find(dat == 253)))=0;
dat(ind2sub(size(dat), find(dat == 254)))=0;
dat(ind2sub(size(dat), find(dat == 255)))=0;
Msk = mask_unprocessed;
Msk.img = dat;
save_untouch_nii(Msk, fullfile(PRD, 'connectivity', 'aparcaseg_2_diff_cortical_only.nii')); 
Msk= fullfile(PRD, 'connectivity', 'aparcaseg_2_diff_cortical_only.nii'); 


%Value of K in 2^K
%K determines the number of nodes.
%Number of nodes is given by 2^K 
% K:                1    2    3    4    5    6    7    8    9    10    11   12
% Number of Nodes:  2    4    8    16   32   64   128  256  512  1024  2048 4096
% K=9;


%DO NOT MODIFY AFTER THIS LINE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Read mask
[hdr,data]=read(Msk);

% check that the whole cortex is in the pcture, otherwise add a o borderline
if sum(sum(data(1, :, :)))~=0
    'zeroed first x'
    data(1, :, :) = zeros(size(data(1, :, :)));
end
if sum(sum(data(end, :, :)))~=0
    'zeroed end x'
    data(end, :, :) = zeros(size(data(end, :, :)));
end
if sum(sum(data(:, 1, :)))~=0
    'zeroed first y'
    data(:, 1, :) = zeros(size(data(:, 1, :)));
end
if sum(sum(data(:, end, :)))~=0
    'zeroed end y'
    data(:, end, :) = zeros(size(data(:, end, :)));
end
if sum(sum(data(:, :,  1)))~=0
    'zeroed first z'
    data(:, :, 1) = zeros(size(data(:, :, 1)));
end
if sum(sum(data(:, :, end)))~=0
    'zeroed end z'
    data(:, :, end) = zeros(size(data(:, :, end)));
end



%Uncomment this to constrain mask to one hemishere:
regions = unique(data);

vol = zeros(size(data));
number_nodes = 0;

for ind_region=2:size(regions, 1) 
     data_curr=zeros(size(data));
    data_curr(ind2sub(size(data), find(data==regions(ind_region))))=1;
    %  current region
    fprintf('Starting recursion for region %d \n', ind_region);
    k=0;
    cell_ind={};
    [cell_ind]=recursive_split(data_curr,cell_ind,k,K);
    vol_curr=zeros(size(data_curr));
    index=randperm(length(cell_ind));
    for i=1:length(cell_ind)
        vol_curr(cell_ind{i})=index(i)+ curr_K*(ind_region-2);
        sz(i)=length(cell_ind{i}); 
    end
    fprintf('Number nodes for region %d: %d, Min size node: %d, Max size node: %d\n',ind_region, length(cell_ind),min(sz),max(sz));
    number_nodes = number_nodes + length(cell_ind);
    vol = vol + vol_curr;
end

% Re add subcortical regions
last_reg = max(vol(:));
mask_unprocessed = load_untouch_nii(fullfile(PRD, '/connectivity/aparcaseg_2_diff.nii.gz')); 
dat = mask_unprocessed.img;
vol(ind2sub(size(dat), find(dat ==  16)))=last_reg +1;
vol(ind2sub(size(dat), find(dat ==   8)))=last_reg +2;
vol(ind2sub(size(dat), find(dat ==  10)))=last_reg +3;
vol(ind2sub(size(dat), find(dat ==  11)))=last_reg +4;
vol(ind2sub(size(dat), find(dat ==  12)))=last_reg +5;
vol(ind2sub(size(dat), find(dat ==  13)))=last_reg +6;
vol(ind2sub(size(dat), find(dat ==  17)))=last_reg +7;
vol(ind2sub(size(dat), find(dat ==  18)))=last_reg +8;
vol(ind2sub(size(dat), find(dat ==  26)))=last_reg +9;
vol(ind2sub(size(dat), find(dat ==  47)))=last_reg +10;
vol(ind2sub(size(dat), find(dat ==  49)))=last_reg +11;
vol(ind2sub(size(dat), find(dat ==  50)))=last_reg +12;
vol(ind2sub(size(dat), find(dat ==  51)))=last_reg +13;
vol(ind2sub(size(dat), find(dat ==  52)))=last_reg +14;
vol(ind2sub(size(dat), find(dat ==  53)))=last_reg +15;
vol(ind2sub(size(dat), find(dat ==  54)))=last_reg +16;
vol(ind2sub(size(dat), find(dat ==  58)))=last_reg +17;

%Output nii file
Out=fullfile(PRD, 'connectivity', ['aparcaseg_2_diff_', num2str(curr_K), '.nii']);


fprintf('Total nodes: %d',number_nodes);
mat2nii(vol,Out,size(data),32,Msk);

%%%%%%%%%%%%%%%%%

% compute centers and orientations
fid = fopen('name_regions.txt');
name_region = textscan(fid, '%s');
fclose(fid);

load('cortical.txt');
cortical(find(cortical==0))=false;
name_region_subcortical = name_region{1}(find(cortical==0));
name_region_subcortical = name_region_subcortical(2:end);
name_region_cortical = name_region{1}(find(cortical==1));

load([PRD,'/', SUBJ_ID, '/connectivity/average_orientations.txt']);
average_orientations = reshape(average_orientations, [88,3]);
average_orientations_subcortical = average_orientations(find(cortical==0), :);
average_orientations_subcortical = average_orientations_subcortical(2:end, :);
average_orientations_cortical = average_orientations(find(cortical==1), :);

name_region = [name_region_cortical; name_region_subcortical];
average_orientations = [average_orientations_cortical; average_orientations_subcortical];
list_region = unique(vol);
list_region = list_region(2:end);
centres = zeros(size(list_region, 1), 4);  
orientation_divided = zeros(size(list_region, 1), 3);

for ind_j=1:(size(list_region, 1)-size(name_region_subcortical, 1))
    list_region(ind_j); 
    [a, b, c] = ind2sub(size(vol), find(vol==list_region(ind_j))); 
    centres(ind_j, 2:4) = [mean(a), mean(b), mean(c)];
    centres(ind_j, 1) = list_region(ceil(ind_j/curr_K)); 
    orientation_divided(ind_j, :) = average_orientations(list_region(ceil(ind_j/curr_K)), :); 
end

for ind_j=(size(list_region, 1)-size(name_region_subcortical, 1)+1):size(list_region, 1) 
    [a, b, c] = ind2sub(size(vol), find(vol==list_region(ind_j))); 
    centres(ind_j, 2:4) = [mean(a), mean(b), mean(c)];
    centres(ind_j, 1) = list_region(ind_j)-(curr_K-1)*70; 
    orientation_divided(ind_j, :) = average_orientations(list_region(ind_j)-(curr_K-1)*70,:); 
end

fid = fopen([PRD, '/', SUBJ_ID, '/connectivity_', num2str(curr_K),'/centres.txt'], 'w'); 
fid2 = fopen([PRD, '/', SUBJ_ID, '/connectivity_', num2str(curr_K),'/average_orientations.txt'], 'w'); 
for i=1:size(list_region, 1)
fprintf(fid, '%s %.3f %.3f %.3f\n', name_region{centres(i, 1)}, centres(i,2:4)'); 
fprintf(fid2, '%.3f %.3f %.3f\n', orientation_divided(i,:)'); 
end
fclose(fid);
fclose(fid2);

% save corr_mat
fid = fopen([PRD, '/connectivity/corr_mat_', num2str(curr_K),'.txt'], 'w'); 
for i=1:size(list_region, 1)
% fprintf(fid, '%d %s\n', list_region(i), name_region{centres(i, 1)}); 
fprintf(fid, '%d %d\n', list_region(i), list_region(i)); 
end
fclose(fid);



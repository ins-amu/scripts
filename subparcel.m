if (~isdeployed) 
addpath('read_and_write_func')   
addpath('uniform_parcellate')
end
curr_K = getenv('curr_K')
K = log(str2num(curr_K))/log(2) 
PRD = getenv('PRD')
SUBJ_ID = getenv('SUBJ_ID')

% get rid of subcortical regions and white matter
mask_unprocessed = load_untouch_nii([PRD, '/connectivity/aparcaseg_2_diff.nii.gz']); 
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
save_untouch_nii(Msk, [PRD, '/connectivity/aparcaseg_2_diff_cortical_only.nii']); 
Msk= [PRD, '/connectivity/aparcaseg_2_diff_cortical_only.nii']; 

%Output nii file
Out=[PRD, '/connectivity/aparcaseg_2_diff_', curr_K, '.nii'];

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

%Uncomment this to constrain mask to one hemishere:

data_1=zeros(size(data));
data_1(ind2sub(size(data), find(data > 1999.9 )))=1;
data_2=zeros(size(data));
data_curr = data;
data_curr(ind2sub(size(data), find(999.9 > data )))=10000;
data_2(ind2sub(size(data), find(data_curr < 1999.9 )))=1;

%Simulated data:
%data=zeros(100,100,5);
%data([5:95],[5:95],2)=1;

% first hemisphere
fprintf('Starting recursion first hemisphere...\n');
fprintf('This may take several minutes\n');
fprintf('(E.g. about 20 mins. to parcellate one hemisphere for K=9)\n');
k=0;
cell_ind={};
[cell_ind]=recursive_split(data_1,cell_ind,k,K-1);

vol_1=zeros(size(data_1));
index=randperm(length(cell_ind));
for i=1:length(cell_ind)
    vol_1(cell_ind{i})=index(i); 
    sz(i)=length(cell_ind{i}); 
end

fprintf('Total nodes: %d, Min size node: %d, Max size node: %d\n',length(cell_ind),min(sz),max(sz));

% second hemisphere
fprintf('Starting recursion second hemisphere...\n');
fprintf('This may take several minutes\n');
fprintf('(E.g. about 20 mins. to parcellate one hemisphere for K=9)\n');
k=0;
cell_ind={};
[cell_ind]=recursive_split(data_2,cell_ind,k,K-1);

vol_2=zeros(size(data_2));
index=randperm(length(cell_ind));
for i=1:length(cell_ind)
    vol_2(cell_ind{i})=index(i)+max(unique(vol_1)); 
    sz(i)=length(cell_ind{i}); 
end

fprintf('Total nodes: %d, Min size node: %d, Max size node: %d\n',length(cell_ind),min(sz),max(sz));

vol = vol_1 + vol_2;
mat2nii(vol,Out,size(data),32,Msk);

%%%%%%%%%%%%%%%%%

% compute centers
list_region = unique(vol);

list_region = list_region(2:end);
centres = zeros(size(list_region, 1), 4);  
for j=1:size(list_region, 1) 
    list_region(j); 
    [a, b, c] = ind2sub(size(vol), find(vol==list_region(j))); 
    centres(j, 2:4) = [mean(a), mean(b), mean(c)];
    centres(j, 1) = int32(list_region(j)); 
end

fid = fopen([PRD, '/', SUBJ_ID, '/connectivity_', num2str(curr_K),'/centres.txt'], 'w'); 
fprintf(fid, '%d %.3f %.3f %.3f\n', centres'); 
fclose(fid);

% save corr_mat
corr_mat = repmat(list_region,1,2);

fid = fopen([PRD, '/connectivity/corr_mat_', curr_K,'.txt'], 'w'); 
fprintf(fid, '%d %d\n', corr_mat'); 
fclose(fid);

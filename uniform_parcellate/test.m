[hdr,data]=read('parcellated.nii');


for i=1:length(unique(data))-1
    sz(i)=length(find(i==data));
end
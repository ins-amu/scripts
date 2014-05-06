dim = 3
curr = 27
figure();
image(squeeze(data(:,:,curr))');


hold on;
for i=5500:8500
    a = tracks.data{i};
    c = -affin(1:3,4);
    d= repmat(c',size(a,1),1);
    e = abs(r*(a+d)')+1;
for j=1:size(e, 1)
if curr-1<e(dim,j) && e(dim,j)<curr+1
scatter(e(1,j),e(2,j));
scatter(e(1,end), e(2,end));
e(:,j);
end
end
end


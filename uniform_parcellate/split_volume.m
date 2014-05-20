function [vol1,vol2]=split_volume(vol)
%Split volume into two contiguous volumes of equal size

ind=find(vol);
sz=size(vol);
ind_surf=get_surface_voxels(vol); 
coor_surf=zeros(length(ind_surf),3);
[coor_surf(:,1),coor_surf(:,2),coor_surf(:,3)]=ind2sub(sz,ind_surf);

%Determine maximally separated seed points on volume surface
%Use show progress if the surface is big
curr_max=0;
if length(ind_surf)<5000
    for i=1:length(ind_surf)
        [max_val,ind_max]=max(((repmat(coor_surf(i,:),length(ind_surf),1)-coor_surf).^2)*[1;1;1]);
        if max_val(1)>curr_max
            curr_max=max_val(1);
            pair=[ind_surf(i),ind_surf(ind_max(1))];
        end
    end
else
    frst=0;
    %Subsample surface voxels
    rp=randperm(length(ind_surf));
    ind_surf=ind_surf(rp(1:5000)); 
    coor_surf_new(:,1)=coor_surf(rp(1:5000),1);
    coor_surf_new(:,2)=coor_surf(rp(1:5000),2);
    coor_surf_new(:,3)=coor_surf(rp(1:5000),3);
    for i=1:length(ind_surf)
        [max_val,ind_max]=max(((repmat(coor_surf_new(i,:),length(ind_surf),1)-coor_surf_new).^2)*[1;1;1]);
        if max_val(1)>curr_max
            curr_max=max_val(1);
            pair=[ind_surf(i),ind_surf(ind_max(1))];
        end
        show_progress(i,length(ind_surf),frst);
        frst=1;
    end 
end
clear coor_surf

%Seed nodes
ind1=find(pair(1)==ind);
ind2=find(pair(2)==ind);

coor=zeros(length(ind),3);
[coor(:,1),coor(:,2),coor(:,3)]=ind2sub(sz,ind);

%Calculate distance between both seed points and all other voxels in volume
dist=zeros(length(ind),2); 
dist(:,1)=((repmat(coor(ind1,:),length(ind),1)-coor).^2)*[1;1;1];
dist(:,2)=((repmat(coor(ind2,:),length(ind),1)-coor).^2)*[1;1;1];

[srt_vals1,srt_ind1]=sort(dist(:,1));
[srt_vals2,srt_ind2]=sort(dist(:,2));
clear dist



srt_ind=[srt_ind1,srt_ind2]';
type=[ones(length(ind),1),ones(length(ind),1)*2]';
type=type(:);
srt_ind=srt_ind(:);

vol1=zeros(sz);
vol2=zeros(sz);
done=zeros(length(srt_ind1),1);
frst=0;
pos1=1;
pos2=1;
stop1=0;
stop2=0;
while ~all(done)

    if stop1==0
    while done(srt_ind1(pos1))
        pos1=pos1+1;
        if pos1>length(srt_ind1)
            stop1=1;
            break
        end
    end
    if pos1<=length(srt_ind1)
        done(srt_ind1(pos1))=1;
        vol1(ind(srt_ind1(pos1)))=1;
    end
    end

    if stop2==0
    while done(srt_ind2(pos2))
        pos2=pos2+1;
        if pos2>length(srt_ind2)
            stop2=1;
            break
        end
    end
    if pos2<=length(srt_ind2)
        done(srt_ind2(pos2))=1;
        vol2(ind(srt_ind2(pos2)))=1;
    end
    end
end

% for i=1:length(type)
%     if ~done(srt_ind(i))
%         done(srt_ind(i))=1;
%         if type(i)==1
%             vol1(ind(srt_ind(i)))=1;
%         else
%             vol2(ind(srt_ind(i)))=1;
%         end
%     end
%     %show_progress(i,length(type),frst);
%     frst=1;
% end

%Assign voxels in volume to one of to seed points, alternating between
%assignments to avoid any bias to one seed point.
% vol1=zeros(sz);
% vol2=zeros(sz);
% flip=1;
% while any(~isinf(srt_ind1))
%     if flip==1
%         while isinf(srt_ind1(1))
%             srt_ind1=srt_ind1(2:length(srt_ind1));
%         end
%         vol1(ind(srt_ind1(1)))=1; 
%         srt_ind2(find(srt_ind1(1)==srt_ind2))=inf; 
%         srt_ind1=srt_ind1(2:length(srt_ind1));
%         flip=0;
%     else
%         while isinf(srt_ind2(1))
%             srt_ind2=srt_ind2(2:length(srt_ind2));
%         end
%         vol2(ind(srt_ind2(1)))=1; 
%         srt_ind1(find(srt_ind2(1)==srt_ind1))=inf; 
%         srt_ind2=srt_ind2(2:length(srt_ind2));
%         flip=1; 
%     end
% end
addpath('iso2mesh/')
PRD = getenv('PRD')
FS = getenv('FS')    

node = load([PRD '/', 'surface/lh_vertices_high.txt']);
d = load([PRD, '/', 'surface/lh_triangles_high.txt']);
face = d + 1;

% we need to ato assign the values to the triangles
% face_size = size(d);
% tri_val = zeros(face_size(1), 1);
% for i=1:face_size(1)
%     if L(face(i,1))==L(face(i,2))
%         tri_val(i) = L(face(i, 1));		
%     elseif L(face(i,1))==L(face(i,3))
%         tri_val(i) = L(face(i, 1));
%     elseif L(face(i,2))==L(face(i,3))
%         tri_val(i) = L(face(i, 2));
%     else
%         tri_val(i) = L(face(i, 1));
%     end
% end

[newno, newfc] = remeshsurf(node, face, 2);
save([PRD, '/', 'surface/lh_vertices_high.txt', 'newno', '-ascii')
f =  newfc(:, 1:3)
save([PRD, '/', 'surface/lh_triangles_high.txt', 'f', '-ascii')
]
% node_newsize = size(newno);
% reg_map = zeros(node_newsize(1),1);
% 
% for  i = 1:node_newsize(1)
% 
%     [row, col] = find(newfc==i);
% 
%     [valuni, induni] = count_unique(newfc(row, 4));
% 
%     [bull, indx] = max(induni);
%     reg_map(i) = bull(indx);
%     
% end

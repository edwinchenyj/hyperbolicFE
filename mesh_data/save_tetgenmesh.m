clear all
close all

% save tetgen mesh from read_tetgenmesh
name = 'scoarsebunny.1';
fs = filesep;
meshname = sprintf('mesh_data%c%s',fs,name);
% if exist([meshname '.mat'], 'file') ~= 2
    [p,~,t] = read_tetgenmesh(meshname);
    if ~isempty(find(t==0))
        t = t' + 1;
    else
        t = t';
    end
    p = p'; p = p(:,1:3);
    temp2 = p(:,2);
    temp3 = p(:,3);
    p(:,3) = temp2;
    p(:,2) = temp3;
    
    [p,t] = fixmesh(p,t);
    nodeM = p./120;
    elem = t;
% else
%     load([meshname '.mat'], 'nodeM', 'elem');
% end
[nodeM,elem] = fixmesh(nodeM,elem); % fix the orientation

%%

for i = 1:size(elem,1)
    T = [nodeM(elem(i,:),:), ones(4,1)];
    assert(det(T)<0)
end


angleX = 145 * pi/180;
angleY = 7 * pi/180;

rxM = [1 0 0; 0 cos(angleX) sin(angleX); 0 -sin(angleX) cos(angleX)];
ryM = [cos(angleY) 0 -sin(angleY); 0 1 0; sin(angleY) 0 cos(angleY)];


for i = 1:size(nodeM,1)
    nodeM(i,:) = nodeM(i,:) * rxM;
    nodeM(i,:) = nodeM(i,:) * ryM;
    
end
simpplot(nodeM,elem,[],[],[],[],[],0.2)


axis equal
grid on
save([meshname '.mat'],'nodeM','elem')
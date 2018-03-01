function test_mass()
clear all
close all


    nodeM = [                % list of xy "node" coordinates
    1, 0
    1, 1
    0, 1
    0, 0                
        ] ;
    
    edge = [                % list of "edges" between nodes
        1, 2                % outer square 
        2, 3
        3, 4
        4, 1
        ] ;

%    [nodeM,etri, ...
%     elem,tnum] = refine2(nodeM,edge) ;
   
hfun = +.1 ;     
      [nodeM,etri, ...
    elem,tnum] = refine2(nodeM,edge,[],[],hfun) ;

elem(:,[1 3]) = elem(:,[3 1]);

trimesh(elem,nodeM(:,1),nodeM(:,2));
N = size(nodeM,1);
M = size(elem,1);
% for i = 1:size(elem,1)
%     T = [nodeM(elem(i,:),:), ones(3,1)];
%     assert(det(T)<0)
% end
obj = staticTriMesh(nodeM, elem);

Y = 100; % Young's modululs
P = 0.4; % Poisson ratio
rho = 1; % density

obj.SetMaterial( Y, P, rho, 1:size(elem,1), 1); % set the tets to be neohookean

obj.M

sum(spdiags(obj.M))
end
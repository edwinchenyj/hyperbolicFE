function initGhostPoints(obj)


% maxA_list = [1 0.1 0.05 0.01]';
maxA = 0.1;

fs = filesep;

type = 'rightangletriangle';

filename = sprintf('simData%ctrimeshes%cmatlabData%c%s maxA%.e barycentric', fs, fs, fs, type, maxA);

load([filename '.mat'], 'Bary','sub_elem_connectivity','num_sub_points','ind_fix','ind_apex','ind_apex2','ind_apex3');
T = 1:3;
    
    
obj.sub_objs = elasticTriObj.empty(obj.NT,0);

for t = 1:obj.NT
    elem = obj.elem(t,:);
    nodeM = obj.nodeM(elem,:);
    
    TR = triangulation(T,nodeM(:,1),nodeM(:,2));
    sub_nodes = barycentricToCartesian(TR,ones(num_sub_points,1),Bary);
    
    obj.sub_objs(t) = elasticTriObj(sub_nodes,sub_elem_connectivity);
    
    obj.sub_objs(t).SetMaterial( obj.Y, obj.P, obj.Rho, 1:size(sub_elem_connectivity,1), obj.material_type);
    obj.sub_objs(t).ind_fix = ind_fix;
    obj.sub_objs(t).ind_apex = ind_apex;
    obj.sub_objs(t).ind_apex2 = ind_apex2;
    obj.sub_objs(t).ind_apex3 = ind_apex3;
    
    obj.sub_objs(t).SetCurrentState(zeros(size(obj.sub_objs(t).X)));
end
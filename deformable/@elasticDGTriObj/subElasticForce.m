function f = subElasticForce(obj)
% compute the force at the tips after relaxing the internal sub points

f = zeros(2*obj.N,1);

maxA = 0.1;


fs = filesep;

type = 'rightangletriangle';


filename = sprintf('simData%ctrimeshes%cmatlabData%c%s maxA%.e barycentric', fs, fs, fs, type, maxA);

load([filename '.mat'], 'Bary','sub_elem_connectivity','num_sub_points','ind_fix','ind_apex','ind_apex2','ind_apex3');
T = [1 2 3];

% prelocate the space for defo
defo = zeros(size(obj.sub_objs(1).X));

options = optimoptions('fsolve','TolFun',1.e-9,'TolX',1.e-9,'Display','final');

for t = 1:obj.NT
    
    i = obj.elem(t, 1); j = obj.elem(t, 2); k = obj.elem(t, 3);
    
    f_new = zeros(2*obj.N,1);
    local_obj = obj.sub_objs(t);
    P = [obj.x(2*(i-1)+1:2*i),...
    obj.x(2*(j-1)+1:2*j),...
    obj.x(2*(k-1)+1:2*k)]';
%     P = [local_obj.x(ind_apex3)';...
%         local_obj.x(ind_apex)';...
%         local_obj.x(ind_apex2)'];
    %     transpose(reshape(local_obj.x(ind_fix),2,3));
    P_original = [obj.X(2*(i-1)+1:2*i),...
    obj.X(2*(j-1)+1:2*j),...
    obj.X(2*(k-1)+1:2*k)]';

    TR = triangulation(T,P_original);
    defo_warm_start = zeros(size(defo));
%     defo_warm_start(ind_fix) = local_obj.x(ind_fix) - local_obj.X(ind_fix);
    defo_warm_start(ind_apex3) = P(1,:)' - local_obj.X(ind_apex3);
    defo_warm_start(ind_apex) = P(2,:)' - local_obj.X(ind_apex);
    defo_warm_start(ind_apex2) = P(3,:)' - local_obj.X(ind_apex2);
    
    
    for i_n = 1:local_obj.N
        PC = local_obj.X(2*(i_n-1)+1:2*i_n)';
        B = cartesianToBarycentric(TR,1,PC);
        defo_warm_start(2*(i_n-1)+1:2*i_n) = B(1) * P(1,:)'...
            + B(2) * P(2,:)' + B(3) * P(3,:)' - PC';
    end
    local_obj.SetCurrentState(defo_warm_start);
    free_defo = defo_warm_start(~ind_fix);
    [dx,fval,exitflag] = fsolve(@local_obj.ElasticForceWFixDx, free_defo, options);
    assert(exitflag > 0);
    
    defo(~ind_fix) = dx;
    defo(ind_apex3) = P(1,:)' - local_obj.X(ind_apex3);
    defo(ind_apex) = P(2,:)' - local_obj.X(ind_apex);
    defo(ind_apex2) = P(3,:)' - local_obj.X(ind_apex2);
    local_obj.SetCurrentState(defo)
    
    force = local_obj.ElasticForce;
    %     tip_force = force(local_obj.ind_fix);
    i = obj.elem(t, 1); j = obj.elem(t, 2); k = obj.elem(t, 3);
    
    f_new(2*(i-1)+1:2*i) = f_new(2*(i-1)+1:2*i)+force(ind_apex3);
    f_new(2*(j-1)+1:2*j) = f_new(2*(j-1)+1:2*j)+force(ind_apex);
    f_new(2*(k-1)+1:2*k) = f_new(2*(k-1)+1:2*k)+force(ind_apex2);
    f = f + f_new;
    
end

end
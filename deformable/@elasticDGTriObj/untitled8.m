function f = subElasticFroce(obj)
% compute the force at the tips after relaxing the internal sub points

f = zeros(2*obj.N,1);


fs = filesep;

type = 'rightangletriangle';

filename = sprintf('simData%ctrimeshes%cmatlabData%c%s maxA%.e barycentric', fs, fs, fs, type, maxA);

load([filename '.mat'], 'Bary','sub_elem_connectivity','num_sub_points','ind_fix');
T = 1:3;

% prelocate the space for Dx
Dx = zeros(size(obj.sub_objs(1).X));

for t = 1:obj.NT
    
    defo_warm_start = zeros(size(Dx));
    free_defo = defo_warm_start(~ind_fix);
    dx = fsolve(@nonlinear_f, free_defo, options);
    
    defo(~ind_fix) = dx;
    obj.SetCurrentState(defo)
    
    tf = 0;
    
end

end
function P = Stress(obj, t)
% calculate the first Piola-Kirchhoff stress tensor of the given tri element
%       t = index of the tri element
material_type = obj.material_type;
tF = obj.F(2*(t-1)+1:2*t,:);
tFINV = obj.FINV(2*(t-1)+1:2*t,:);

mu = obj.mu;
lambda = obj.lambda;
if material_type == 1 % neo-hookean
    J = det(tF);
    P = mu *(tF - tFINV') + lambda * log(J) * tFINV;
elseif material_type == 2 % linear elasticity
    P = mu*(tF + tF' - 2*obj.Iv) + lambda*trace(tF - obj.Iv)*obj.Iv;
else
    error('Unexpect error. No material type specified')
end
end
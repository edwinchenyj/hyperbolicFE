function P = Stress(obj, t)
% calculate the first Piola-Kirchhoff stress tensor of the given tri element
%       t = index of the tri element
material_type = obj.elemMaterialType(t);
tF = obj.F(2*(t-1)+1:2*t,:);
tFINV = obj.FINV(2*(t-1)+1:2*t,:);

mu = obj.mu(t);
lambda = obj.lambda(t);
if material_type == 1 % neo-hookean
    J = det(tF);
    P = mu *(tF - tFINV') + lambda * log(J) * tFINV;
elseif material_type == 2 % linear elasticity
    P = mu*(tF + tF' - 2*obj.I2) + lambda*trace(tF - obj.I2)*obj.I2;
else
    error('Unexpect error. No material type specified')
end
end
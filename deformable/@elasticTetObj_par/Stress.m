function P = Stress(obj, t)
% calculate the first Piola-Kirchhoff stress tensor of the given tet element
%       t = index of the tetrahedron element
type = obj.elemMaterialType(t);
tF = obj.F(3*(t-1)+1:3*t,:);
tFINV = obj.FINV(3*(t-1)+1:3*t,:);

mu = obj.mu(t);
lambda = obj.lambda(t);
if type == 1 % neo-hookean
    J = det(tF);
    P = mu *(tF - tFINV') + lambda * log(J) * tFINV;
elseif type == 2 % linear elasticity
    P = mu*(tF + tF' - 2*obj.I3) + lambda*trace(tF - obj.I3)*obj.I3;
else
    error('Unexpect error. No material type specified')
end
end
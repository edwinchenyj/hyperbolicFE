function dP = StressDifferential(obj, t, dF)
% calculate the stress differential of the given tet element
% and pre-calculated differential of the deformation gradient
type = obj.elemMaterialType(t);
tF = obj.F(2*(t-1)+1:2*t,:);
tFINV = obj.FINV(2*(t-1)+1:2*t,:);
dP = zeros(2,2);

mu = obj.mu(t);
lambda = obj.lambda(t);
if type == 1 % neo-hookean
    dP = mu * dF+ mu * (tFINV * dF * tFINV)' - lambda * log(det(tF)) * (tFINV * dF * tFINV)'...
                + lambda * trace(tFINV * dF) * (tFINV)';
%     disp('dP= [')
%     disp(dP)
%     disp('];')
elseif type == 2 % linear elasticity
    dP = mu * (dF + dF') + lambda * trace(dF) * obj.I2;
else
    error('Unexpect error. No material type specified')
end
end
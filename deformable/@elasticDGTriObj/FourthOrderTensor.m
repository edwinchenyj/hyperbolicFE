function C = FourthOrderTensor(obj, t)
% calculate the forth order tensor corresponding to the first
% Piola-Kirchhoff stress tensor of the given tri element 
% C = \frac{\partial^2}{\partial x^T \partial x}E
%       t = index of the tri element
mu = obj.mu(t);
lambda = obj.lambda(t);
tF = obj.F(2*(t-1)+1:2*t,:);
tFINV = obj.FINV(2*(t-1)+1:2*t,:);

% the fourth order tensor
if (obj.elemMaterialType(t) == 1)
    % for Neo-hookean
    C = mu * obj.I4 + mu * obj.K44* kron(tFINV',tFINV)...
        - lambda * (log(det(tF))*obj.K44*kron(tFINV',tFINV))...
        + lambda*(obj.K44*(tFINV(:)*reshape(transpose(tFINV),1,4)));
elseif (obj.elemMaterialType(t) == 2)
    % for linear elasticity
    C = mu * (obj.I4 + obj.K44) + lambda * (obj.K44 * obj.I2(:)*obj.I2(:)');
elseif (obj.elemMaterialType(t) == 3)
    
end
end
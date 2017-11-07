function C = FourthOrderTensor(obj, t)
% calculate the forth order tensor corresponding to the first
% Piola-Kirchhoff stress tensor of the given tri element 
% C = \frac{\partial^2}{\partial x^T \partial x}E
%       t = index of the tri element
mu = obj.mu;
lambda = obj.lambda;
tF = obj.F(2*(t-1)+1:2*t,:);
tFINV = obj.FINV(2*(t-1)+1:2*t,:);

% the fourth order tensor
if (obj.material_type == 1)
    % for Neo-hookean
    C = mu * obj.Im + mu * obj.Kmm* kron(tFINV',tFINV)...
        - lambda * (log(det(tF))*obj.Kmm*kron(tFINV',tFINV))...
        + lambda*(obj.Kmm*(tFINV(:)*reshape(transpose(tFINV),1,4)));
elseif (obj.material_type == 2)
    % for linear elasticity
    C = mu * (obj.Im + obj.Kmm) + lambda * (obj.Kmm * obj.Iv(:)*obj.Iv(:)');
elseif (obj.material_type == 3)
    
end
end
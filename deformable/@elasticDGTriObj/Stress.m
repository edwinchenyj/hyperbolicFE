function P = Stress(obj, t)
% calculate the first Piola-Kirchhoff stress tensor of the given tri element
%       t = index of the tri elmu = obj.mu(t);
            lambda = obj.lambda(t);
            tT = obj.T(4*(t-1)+1:4*t,:);
            W = obj.W(t);
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
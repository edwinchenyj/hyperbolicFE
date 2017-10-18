function f = ElasticForce(obj)
% compute the elastic force under the current deformation (3N by 1)
% can be calculated with current state as the sole input...

f = zeros(2*obj.N,1);
switch obj.elemMaterialType(1) % TODO: change the material type initialization
    case 1
        
        for t = 1:obj.NT
            
            f_new = zeros(2*obj.N,1);
            
            type = obj.elemMaterialType(t);
            tF = obj.F(2*(t-1)+1:2*t,:);
            tFINV = obj.FINV(2*(t-1)+1:2*t,:);
            
            mu = obj.mu(t);
            lambda = obj.lambda(t);
            
            
            
            
            if type == 1 % neo-hookean
                J = det(tF);
                P = mu *(tF - tFINV') + lambda * log(J) * tFINV';
            elseif type == 2 % linear elasticity
                P = mu*(tF + tF' - 2*obj.I2) + lambda*trace(tF - obj.I2)*obj.I2;
            elseif type == 3
                E = 1/2 * (F'*F - I);
                P = F*(2*mu * E + lambda * trace(E) * I);
            else
                error('Unexpect error. No material type specified')
            end
            
            H = -obj.W(t) * P * (obj.DmINV(2*(t-1)+1:2*t,:)');
            i = obj.elem(t, 1); j = obj.elem(t, 2); k = obj.elem(t, 3);
            
            f_new(2*(i-1)+1:2*i) = f_new(2*(i-1)+1:2*i)+H(:,1);
            f_new(2*(j-1)+1:2*j) = f_new(2*(j-1)+1:2*j)+H(:,2);
            f_new(2*(k-1)+1:2*k) = f_new(2*(k-1)+1:2*k) - H(:,1) - H(:,2);
            f = f + f_new;
        end
    case 2
        if isempty(obj.K0)
            
            for t = 1:obj.NT
                
                mu = obj.mu(t);
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
                %     obj.C = C;
                %     disp('dFdx = [')
                %     disp(tT)
                %     disp('];')
                %     obj.dFdx = tT;
                %     disp('C= [')
                %     disp(C)
                %     disp('];')
                
                
                % element stiffness matrix
                
                Kt = W * tT'*C*tT;
                Kt = 1/2 * (Kt + Kt');
                
                for ti = 1:3
                    for tj = 1:3
                        %            sA(index:index+3) = Kt(obj.IndK(3*(ti-1) + tj,:));
                        
                        % same as:
                        sA(index:index+3) = reshape(Kt(2*(ti-1)+1:2*ti,2*(tj-1)+1:2*tj),4,1);
                        index = index + 4;
                    end
                end
            end
            
            % global stiffness matrix
            obj.K0 = sparse(obj.ii, obj.jj, sA, 2*obj.N, 2*obj.N);
            
        end
        f = -obj.K0 * (obj.x - obj.X);
end
end
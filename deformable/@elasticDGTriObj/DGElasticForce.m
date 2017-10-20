function DGf_k = DGElasticForce(obj)
% compute the elastic force under the current deformation (3N by 1)
% can be calculated with current state as the sole input...

DGf_k = zeros(2*obj.DGN,1);

switch obj.elemMaterialType(1)
    case 1
        for t = 1:obj.NT
            
            DGf_k_new = zeros(2*obj.DGN,1);
            
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
            i = obj.DGelem(t, 1); j = obj.DGelem(t, 2); k = obj.DGelem(t, 3);
            
            DGf_k_new(2*(i-1)+1:2*i) = DGf_k_new(2*(i-1)+1:2*i)+H(:,1);
            DGf_k_new(2*(j-1)+1:2*j) = DGf_k_new(2*(j-1)+1:2*j)+H(:,2);
            DGf_k_new(2*(k-1)+1:2*k) = DGf_k_new(2*(k-1)+1:2*k) - H(:,1) - H(:,2);
            DGf_k = DGf_k + DGf_k_new;
        end
        
        DGf_k = DGf_k + obj.DGInterfaceElasticForce;
        
    case 2
        if isempty(obj.DGK0)
            index = 1;

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
            DGK_k = sparse(obj.DGElement_ii, obj.DGElement_jj, sA, 2*obj.DGN, 2*obj.DGN);
            
            % DGK_interface = sparse(size(obj.DGX,1),size(obj.DGX,1));
            % traverse the interface
            obj.DGK0 = DGK_k + obj.DGInterfaceStiffnessMatrix;
        end
        DGf_k = -obj.DGK0 * (obj.DGx - obj.DGX);
        
end
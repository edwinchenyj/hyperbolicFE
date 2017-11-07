function DGf_k = ElasticForce(obj)
% compute the elastic force under the current deformation (3N by 1)
% can be calculated with current state as the sole input...

DGf_k = zeros(2*obj.N,1);

switch obj.material_type
    case 1
        for t = 1:obj.NT
            
            DGf_k_new = zeros(2*obj.N,1);
            
            type = obj.material_type;
            tF = obj.F(2*(t-1)+1:2*t,:);
            tFINV = obj.FINV(2*(t-1)+1:2*t,:);
            
            mu = obj.mu(t);
            lambda = obj.lambda(t);
            
            
            
            
            if type == 1 % neo-hookean
                J = det(tF);
                P = mu *(tF - tFINV') + lambda * log(J) * tFINV';
            elseif type == 2 % linear elasticity
                P = mu*(tF + tF' - 2*obj.Iv) + lambda*trace(tF - obj.Iv)*obj.Iv;
            elseif type == 3
                E = 1/2 * (F'*F - I);
                P = F*(2*mu * E + lambda * trace(E) * I);
            else
                error('Unexpect error. No material type specified')
            end
            
            H = -obj.W(t) * P * (obj.DmINV(2*(t-1)+1:2*t,:)');
            i = obj.elem(t, 1); j = obj.elem(t, 2); k = obj.elem(t, 3);
            
            DGf_k_new(2*(i-1)+1:2*i) = DGf_k_new(2*(i-1)+1:2*i)+H(:,1);
            DGf_k_new(2*(j-1)+1:2*j) = DGf_k_new(2*(j-1)+1:2*j)+H(:,2);
            DGf_k_new(2*(k-1)+1:2*k) = DGf_k_new(2*(k-1)+1:2*k) - H(:,1) - H(:,2);
            DGf_k = DGf_k + DGf_k_new;
        end
        
        DGf_k = DGf_k + obj.DGInterfaceElasticForce;
        
    case 2
        if isempty(obj.K0)
            index = 1;

            for t = 1:obj.NT
                
                mu = obj.mu(t);
                lambda = obj.lambda(t);
                tT = obj.T(4*(t-1)+1:4*t,:);
                W = obj.W(t);
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
            DGK_k = sparse(obj.DGElement_ii, obj.DGElement_jj, sA, 2*obj.N, 2*obj.N);
            
            % DGK_interface = sparse(size(obj.DGX,1),size(obj.DGX,1));
            % traverse the interface
            obj.K0 = DGK_k + obj.DGInterfaceStiffnessMatrix;
        end
        DGf_k = -obj.K0 * (obj.x - obj.X);
        
end
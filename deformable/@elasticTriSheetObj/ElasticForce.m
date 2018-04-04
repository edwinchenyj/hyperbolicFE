function f = ElasticForce(obj)
% compute the elastic force under the current deformation (3N by 1)
% can be calculated with current state as the sole input...

f = zeros(2*obj.N,1);
switch obj.material_type
    case 1
        
        for t = 1:obj.NT
            
            f_new = zeros(2*obj.N,1);
            
            type = obj.material_type;
            tF = obj.F(2*(t-1)+1:2*t,:);
            tFINV = obj.FINV(2*(t-1)+1:2*t,:);
            
            mu = obj.mu;
            lambda = obj.lambda;
            
            
            
            
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
            
            f_new(2*(i-1)+1:2*i) = f_new(2*(i-1)+1:2*i)+H(:,1);
            f_new(2*(j-1)+1:2*j) = f_new(2*(j-1)+1:2*j)+H(:,2);
            f_new(2*(k-1)+1:2*k) = f_new(2*(k-1)+1:2*k) - H(:,1) - H(:,2);
            f = f + f_new;
        end
    case 2
        if isempty(obj.K0)
            index = 1;
            for t = 1:obj.NT
                
                mu = obj.mu;
                lambda = obj.lambda;
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
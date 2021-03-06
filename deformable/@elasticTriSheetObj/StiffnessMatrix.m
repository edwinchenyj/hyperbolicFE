function K = StiffnessMatrix(obj)
% construct and return the stiffness matrix under the current deformation,
% this can be calculated using current deformation state as the only input

index = 1;
sA = zeros(size(obj.ii));
switch obj.material_type % TODO: change the material type initialization
    case 1
        
        for t = 1:obj.NT
            
            mu = obj.mu;
            lambda = obj.lambda;
            tT = obj.T(4*(t-1)+1:4*t,:);
            W = obj.W(t);
            tF = obj.F(2*(t-1)+1:2*t,:);
            tFINV = obj.FINV(2*(t-1)+1:2*t,:);
            
            % the fourth order tensor
                % for Neo-hookean
                C = (mu * obj.Im + lambda * reshape(tFINV',[],1)*reshape(tFINV',[],1)' - (lambda * log(det(tF)) - mu) * kron(tFINV,tFINV') * obj.Kmm);
%                 C = mu * obj.Im + mu * obj.Kmm* kron(tFINV',tFINV)...
%                     - lambda * (log(det(tF))*obj.Kmm*kron(tFINV',tFINV))...
%                     + lambda*(obj.Kmm*(tFINV(:)*reshape(transpose(tFINV),1,4)));
%             C = mu * obj.Im + mu * obj.Kmm* kron(tFINV',tFINV)...
%                     - lambda * (log(det(tF))*obj.Kmm*kron(tFINV',tFINV))...
%                     + lambda*(obj.Kmm*(tFINV(:)*reshape(transpose(tFINV),1,4)));
            
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
        K = sparse(obj.ii, obj.jj, sA, 2*obj.N, 2*obj.N);
    case 2
        if isempty(obj.K0)
            
            for t = 1:obj.NT
                
                mu = obj.mu;
                lambda = obj.lambda;
                tT = obj.T(4*(t-1)+1:4*t,:);
                W = obj.W(t);
                tF = obj.F(2*(t-1)+1:2*t,:);
                tFINV = obj.FINV(2*(t-1)+1:2*t,:);
                
                % the fourth order tensor
                    % for linear elasticity
                    C = mu * (obj.Im + obj.Kmm) + lambda * (obj.Kmm * obj.Iv(:)*obj.Iv(:)');
                
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
        K = obj.K0;
        
end
end
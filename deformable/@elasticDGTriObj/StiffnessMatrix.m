function DGK = StiffnessMatrix(obj)
% construct and return the element DG stiffness matrix under the current deformation,
% this can be calculated using current deformation state as the only input

index = 1;
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
            if (obj.material_type == 1)
                % for Neo-hookean
                C = mu * obj.Im + mu * kron(tFINV,tFINV') * obj.Kmm...
                    - lambda * (log(det(tF))*kron(tFINV,tFINV')*obj.Kmm)...
                    + lambda*(reshape(tFINV',[],1) * reshape(tFINV',[],1)');
%                 C2 = mu * obj.I4 + mu * obj.K44* kron(tFINV',tFINV)...
%                     - lambda * (log(det(tF))*obj.K44*kron(tFINV',tFINV))...
%                     + lambda*(obj.K44*(tFINV(:)*reshape(transpose(tFINV),1,4)));

            elseif (obj.material_type == 2)
                % for linear elasticity
                C = mu * (obj.Im + obj.Kmm) + lambda * (obj.Kmm * obj.Iv(:)*obj.Iv(:)');
            elseif (obj.material_type == 3)
                
            end
            
            % element stiffness matrix
            
            Kt = W * tT'*C*tT;
            
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
        
        % DGK_interface = sparse(size(obj.X,1),size(obj.X,1));
        % traverse the interface
        DGK = DGK_k + obj.DGInterfaceStiffnessMatrix;
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
            
            % DGK_interface = sparse(size(obj.X,1),size(obj.X,1));
            % traverse the interface
            obj.K0 = DGK_k + obj.DGInterfaceStiffnessMatrix;
            DGK = obj.K0;
            
        else
            DGK = obj.K0;
        end
end
nonsym = full(max(max(DGK - DGK')));
assert(nonsym < 1e-6,sprintf('stiffness matrix should be symmetric. the error is  %d',nonsym)); 
end
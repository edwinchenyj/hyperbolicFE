function DGK = DGStiffnessMatrix(obj)
% construct and return the element DG stiffness matrix under the current deformation,
% this can be calculated using current deformation state as the only input

index = 1;
switch obj.elemMaterialType(1) % TODO: change the material type initialization
    case 1
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
                C = mu * obj.I4 + mu * kron(tFINV,tFINV') * obj.K44...
                    - lambda * (log(det(tF))*kron(tFINV,tFINV')*obj.K44)...
                    + lambda*(reshape(tFINV',[],1) * reshape(tFINV',[],1)');
%                 C2 = mu * obj.I4 + mu * obj.K44* kron(tFINV',tFINV)...
%                     - lambda * (log(det(tF))*obj.K44*kron(tFINV',tFINV))...
%                     + lambda*(obj.K44*(tFINV(:)*reshape(transpose(tFINV),1,4)));

            elseif (obj.elemMaterialType(t) == 2)
                % for linear elasticity
                C = mu * (obj.I4 + obj.K44) + lambda * (obj.K44 * obj.I2(:)*obj.I2(:)');
            elseif (obj.elemMaterialType(t) == 3)
                
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
        DGK_k = sparse(obj.DGElement_ii, obj.DGElement_jj, sA, 2*obj.DGN, 2*obj.DGN);
        
        % DGK_interface = sparse(size(obj.DGX,1),size(obj.DGX,1));
        % traverse the interface
        DGK = DGK_k + obj.DGInterfaceStiffnessMatrix;
    case 2
        if isempty(obj.DGK0)
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
            DGK_k = sparse(obj.DGElement_ii, obj.DGElement_jj, sA, 2*obj.DGN, 2*obj.DGN);
            
            % DGK_interface = sparse(size(obj.DGX,1),size(obj.DGX,1));
            % traverse the interface
            obj.DGK0 = DGK_k + obj.DGInterfaceStiffnessMatrix;
            DGK = obj.DGK0;
            
        else
            DGK = obj.DGK0;
        end
end
nonsym = full(max(max(DGK - DGK')));
assert(nonsym < 1e-6,sprintf('stiffness matrix should be symmetric. the error is  %d',nonsym)); 
end
function DGK_k = DGElementStiffnessMatrix(obj)
% construct and return the element DG stiffness matrix under the current deformation,
% this can be calculated using current deformation state as the only input

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
end
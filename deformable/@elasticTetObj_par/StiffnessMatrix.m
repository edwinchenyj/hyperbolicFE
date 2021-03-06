function K = StiffnessMatrix(obj)
% construct and return the stiffness matrix under the current deformation,
% this can be calculated using current deformation state as the only input

index = 1;
for t = 1:obj.NT
    
    mu = obj.mu;
    lambda = obj.lambda;
    tT = obj.T(9*(t-1)+1:9*t,:);
    W = obj.W(t);
    tF = obj.F(3*(t-1)+1:3*t,:);
    tFINV = obj.FINV(3*(t-1)+1:3*t,:);
    
    if obj.isCorotatedLinear
        C = mu * (obj.Im + obj.Kmm) + lambda * (obj.Kmm * obj.Iv(:)*obj.Iv(:)');
        Kt = W * tT'*C*tT;
        Kt = 1/2 * (Kt + Kt');
        
        % warped stiffness
        [u,s,v] = svd(tF);
        
        % store elementwise singular vectors
        obj.U(3*(t-1)+1:3*t,:) = u;
        obj.V(3*(t-1)+1:3*t,:) = v;
        obj.S(3*(t-1)+1:3*t) = diag(s);
        
        tR = v'*u;
        % store elementwise rotation
        obj.R(3*(t-1)+1:3*t,:) = tR;
        
        R4 = repmat(tR,1,4);
        i     = repmat(reshape(1:12,3,4),3,1);
        j     = repmat(1:12,3,1);
        R     = sparse(i(:),j(:),R4(:));
        Kt = R*Kt*R';
    else
        % the fourth order tensor
        if (obj.material_type == 1)
            % for Neo-hookean
            C = mu * obj.Im + mu * obj.Kmm* kron(tFINV',tFINV)...
                - lambda * (log(det(tF))*obj.Kmm*kron(tFINV',tFINV))...
                + lambda*(obj.Kmm*(tFINV(:)*reshape(transpose(tFINV),1,9)));
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
    end
    
    for ti = 1:4
        for tj = 1:4
            sA(index:index+8) = Kt(obj.IndK(4*(ti-1) + tj,:));
            
            % same as:
            % sA(index:index+8) = reshape(Kt(3*(ti-1)+1:3*ti,3*(tj-1)+1:3*tj),9,1);
            index = index + 9;
        end
    end
end

% global stiffness matrix
K = sparse(obj.ii, obj.jj, sA, 3*obj.N, 3*obj.N);
end
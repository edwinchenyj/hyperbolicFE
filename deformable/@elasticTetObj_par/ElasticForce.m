function f = ElasticForce(obj)
% compute the elastic force under the current deformation (3N by 1)
% can be calculated with current state as the sole input...

f = zeros(3*obj.N,1);

for t = 1:obj.NT
    
    f_new = zeros(3*obj.N,1);
    
    type = obj.material_type;
    tF = obj.F(3*(t-1)+1:3*t,:);
    tFINV = obj.FINV(3*(t-1)+1:3*t,:);
    
    mu = obj.mu;
    lambda = obj.lambda;
    
    
    if obj.isCorotatedLinear
        tT = obj.T(9*(t-1)+1:9*t,:);
        W = obj.W(t);
        T_nodeM = obj.nodeM(obj.elem(t,:),:);
        T_node = obj.node(obj.elem(t,:),:);
        T_x = T_node';
        T_x = T_x(:); % vectorized nodal position
        
        T_X = T_nodeM';
        T_X = T_X(:); % vectorized undeformed nodal position
        
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
        rI     = repmat(reshape(1:12,3,4),3,1); % row indices
        cI     = repmat(1:12,3,1); % column indices
        R     = sparse(rI(:),cI(:),R4(:));
        i = obj.elem(t, 1); j = obj.elem(t, 2); k = obj.elem(t, 3); l = obj.elem(t, 4);
        
        
        tf = -R*Kt*(R'*T_x-T_X);
        
        f_new(3*(i-1)+1:3*i) = f_new(3*(i-1)+1:3*i)+tf(1:3);
        f_new(3*(j-1)+1:3*j) = f_new(3*(j-1)+1:3*j)+tf(4:6);
        f_new(3*(k-1)+1:3*k) = f_new(3*(k-1)+1:3*k)+tf(7:9);
        f_new(3*(l-1)+1:3*l) = f_new(3*(l-1)+1:3*l)+tf(10:12);
        f = f + f_new;
    else
    

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
    
    H = -obj.W(t) * P * (obj.DmINV(3*(t-1)+1:3*t,:)');
    i = obj.elem(t, 1); j = obj.elem(t, 2); k = obj.elem(t, 3); l = obj.elem(t, 4);
    
    f_new(3*(i-1)+1:3*i) = f_new(3*(i-1)+1:3*i)+H(:,1);
    f_new(3*(j-1)+1:3*j) = f_new(3*(j-1)+1:3*j)+H(:,2);
    f_new(3*(k-1)+1:3*k) = f_new(3*(k-1)+1:3*k)+H(:,3);
    f_new(3*(l-1)+1:3*l) = f_new(3*(l-1)+1:3*l) - H(:,1) - H(:,2) - H(:,3);
    f = f + f_new;
    end

end
end
function DGf_k = DGElementElasticForce(obj)
% compute the elastic force under the current deformation (3N by 1)
% can be calculated with current state as the sole input...

DGf_k = zeros(2*obj.N,1);

for t = 1:obj.NT
    
    DGf_k_new = zeros(2*obj.N,1);
    
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
    
    DGf_k_new(2*(i-1)+1:2*i) = DGf_k_new(2*(i-1)+1:2*i)+H(:,1);
    DGf_k_new(2*(j-1)+1:2*j) = DGf_k_new(2*(j-1)+1:2*j)+H(:,2);
    DGf_k_new(2*(k-1)+1:2*k) = DGf_k_new(2*(k-1)+1:2*k) - H(:,1) - H(:,2);
    DGf_k = DGf_k + DGf_k_new;
end


end
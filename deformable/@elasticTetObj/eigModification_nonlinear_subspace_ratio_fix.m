function [K, f] = eigModification_nonlinear_subspace_ratio_fix(obj)

if (obj.material_type == 2) 
    if isempty(obj.modifiedK)
    K = obj.StiffnessMatrix;
    M = obj.M;
    
    M = M(obj.indLogical,obj.indLogical);
    K = K(obj.indLogical,obj.indLogical);
    [V,D] = eig(full(K),full(M));
    eDiag = diag(D);
    [eDiag, permutation] = sort(eDiag);
    V = V(:,permutation);
    
    for i = 1:length(obj.eig_targets)
        eDiag(i) = obj.eig_targets(i);
    end
    Minv_K_new = V*diag(eDiag)/V;
    obj.modifiedK = M * Minv_K_new;
    end
K = obj.modifiedK;
f = -K * (obj.x - obj.X);

else
    K = obj.StiffnessMatrix;
    M = obj.M;
    
    M = M(obj.indLogical,obj.indLogical);
    K = K(obj.indLogical,obj.indLogical);
    
    [V_s,D] = eigs(K,M,length(obj.eig_targets),'sm');
    [low_eig, permutation_indices] = sort(diag(D));
%     obj.eig_ratios = obj.eig_targets./low_eig;
    V_s = V_s(:,permutation_indices);
    rescaled_K = K + M * V_s * (diag(obj.eig_ratios) - eye(length(obj.eig_ratios))) * (V_s')*K*V_s * (V_s') * M;

    Eforce = obj.ElasticForce;
    Eforce = Eforce(obj.indLogical);
    rescaled_f =  Eforce + M * V_s * (diag(obj.eig_ratios) - eye(length(obj.eig_ratios))) * (V_s') * Eforce;
    
    K = rescaled_K;
    f = rescaled_f;
end
end
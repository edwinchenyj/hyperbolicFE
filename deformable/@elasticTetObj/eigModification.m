function K = eigModification(obj,eigenvalues)

if isempty(obj.modifiedK)
    K = obj.StiffnessMatrix;
    M = obj.M;
    
    M = M(obj.indLogical,obj.indLogical);
    K = K(obj.indLogical,obj.indLogical);
    [V_s,D] = eigs(K,M,length(obj.eig_targets),'sm');
    [low_eig, permutation_indices] = sort(diag(D));
    V_s = V_s(:,permutation_indices);
    obj.eig_ratios = eigenvalues./low_eig;
    rescaled_K = K + M * V_s * (diag(obj.eig_ratios) - eye(length(obj.eig_ratios))) * (V_s')*K*V_s * (V_s') * M;

    obj.modifiedK = rescaled_K;
    
end
K = obj.modifiedK;
end
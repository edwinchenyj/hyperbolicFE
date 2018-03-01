function [K, f] = eigModification_nonlinear_subspace(obj)

if (obj.material_type == 2)
    
    if isempty(obj.modifiedK)
        K = obj.StiffnessMatrix;
        M = obj.M;
        obj.modifiedK = K; % store the K for the original size
        
        M = M(obj.indLogical,obj.indLogical);
        K = K(obj.indLogical,obj.indLogical);
        [V_s,D] = eigs(K,M,length(obj.eig_targets),'sm');
        [low_eig, permutation_indices] = sort(diag(D));
        V_s = V_s(:,permutation_indices);
        obj.eig_ratios = obj.eig_targets./low_eig;

        rescaled_K = K + M * V_s * (diag(obj.eig_ratios) - eye(length(obj.eig_ratios))) * (diag(low_eig)) * (V_s') * M;
%         rescaled_K = K;
        obj.modifiedK(obj.indLogical,obj.indLogical) = rescaled_K;
        
    end
    K = obj.modifiedK(obj.indLogical,obj.indLogical);
    f = -obj.modifiedK * (obj.x - obj.X);
    f = f(obj.indLogical);
else
    K = obj.StiffnessMatrix;
    M = obj.M;
    
    M = M(obj.indLogical,obj.indLogical);
    K = K(obj.indLogical,obj.indLogical);
    
    [V_s,D] = eigs(K,M,length(obj.eig_targets),'sm');
    [low_eig, permutation_indices] = sort(diag(D));
    obj.eig_ratios = obj.eig_targets./low_eig;
    V_s = V_s(:,permutation_indices);
%     rescaled_K = K + M * V_s * (diag(obj.eig_ratios) - eye(length(obj.eig_ratios))) * (V_s')*K*V_s * (V_s') * M;
    rescaled_K = K + M * V_s * (diag(obj.eig_ratios) - eye(length(obj.eig_ratios))) * (diag(low_eig)) * (V_s') * M;
    
    Eforce = obj.ElasticForce;
    Eforce = Eforce(obj.indLogical);
    rescaled_f =  Eforce + M * V_s * (diag(obj.eig_ratios) - eye(length(obj.eig_ratios))) * (V_s') * Eforce;
    
    K = rescaled_K;
    f = rescaled_f;
end
end
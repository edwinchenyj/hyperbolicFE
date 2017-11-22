function energy= ElasticEnergy(obj)

energy  = 0;

elastic_energy = 0;
for t = 1:obj.NT
    
    tF = obj.F(2*(t-1)+1:2*t,:);
    mu = obj.mu;
    lambda = obj.lambda;
    W = obj.W(t);
    
    if (obj.material_type == 1)
        % for Neo-hookean
%         [~,s,~] = svd(tF);
        I1 = trace(tF'*tF);
        I3 = (det(tF))^2;
        
        tEnergy = W * (mu/2 * (I1 - log(I3) - 2) + lambda/4 * ((log(I3))^2) );
        
    elseif (obj.material_type == 2)
        % for linear elasticity
        linearStrain = 1/2 * (tF + tF') - eye(2);
        tEnergy = W * (mu * trace(linearStrain'*linearStrain) + lambda/2 * (trace(linearStrain)^2));
        
    elseif (obj.material_type == 3)
        
    end
    
    elastic_energy = elastic_energy + tEnergy;
    
end

energy = elastic_energy;

end
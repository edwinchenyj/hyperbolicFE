function energy= totalEnergy(obj)

energy  = 0;

elastic_energy = 0;
for t = 1:obj.NT
    
    tF = obj.F(3*(t-1)+1:3*t,:);
    mu = obj.mu(t);
    lambda = obj.lambda(t);
    W = obj.W(t);
    
    if (obj.elemMaterialType(t) == 1)
        % for Neo-hookean
%         [~,s,~] = svd(tF);
        I1 = trace(tF'*tF);
        I3 = (det(tF))^2;
        
        tEnergy = W * (mu/2 * (I1 - log(I3) - 3) + lambda/8 * ((log(I3))^2) );
        
    elseif (obj.elemMaterialType(t) == 2)
        % for linear elasticity
        linearStrain = 1/2 * (tF + tF') - eye(3);
        tEnergy = W * (mu * trace(linearStrain'*linearStrain) + lambda/2 * (trace(linearStrain)^2));
        
    elseif (obj.elemMaterialType(t) == 3)
        
    end
    
    elastic_energy = elastic_energy + tEnergy;
    
end
diagMass = diag(obj.M);
nodeMass = diagMass(1:3:end);
nodeV = reshape(obj.v,3,obj.N);
nodeV2 = sum(nodeV.^2,1);
nodeV2 = nodeV2';
kinetic = 1/2 * nodeMass' * nodeV2;

z = obj.x;
z = z(3:3:end);
potential = 9.8 * nodeMass' * z;

energy = elastic_energy + kinetic + potential;

end
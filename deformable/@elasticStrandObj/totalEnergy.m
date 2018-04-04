function energy= totalEnergy(obj)

energy  = 0;

elastic_energy = 0;
for s = 1:obj.NS
    
    F = obj.F(3*(s-1)+1:3*s,:);    
    ks = obj.ks(s);
    l0 = obj.W(s);
    
    
    elastic_energy = elastic_energy + 1/8 * l0 * ks * (F'*F - 1)^2 * l0;
    
    
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
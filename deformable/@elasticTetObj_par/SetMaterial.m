function SetMaterial(obj,Y, P, Rho, type, a, b)
% Set the material type of some specified elements
%   elem: the vector of elements specified
%   type: 1 for neo-hookean, 2 for linear, and more to come...

assert(type == 1 || type == 2);

obj.a = a;
obj.b = b;

obj.M = sparse(size(obj.M,1),size(obj.M,2));
obj.mu = Y ./ ( 2 * (1 + P) );
obj.lambda = ( Y .* P ) ./ ( (1 + P) .* (1 - 2 * P) );
obj.rho = Rho;
obj.Y = Y;
obj.P = P;
obj.Rho = Rho;
obj.material_type = type;
for t = 1:obj.NT
%     assert(obj.elemMaterialType(elem(t)) == 0); % check the element haven't been initialized to any material type

    for e = obj.elem(t,:)
        for mi = (e-1)*3+1:e*3
            obj.M(mi,mi) = obj.M(mi,mi)+obj.W(t)/4 * Rho;
        end
    end
end

% after setting the material, calculate the rest_stiffness
obj.SetCurrentState(zeros(obj.Dim*obj.N,1));
obj.rest_stiffness = obj.StiffnessMatrix;

end
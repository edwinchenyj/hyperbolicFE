function SetMaterial(obj,Y, P, Rho, mtype, a, b)
% Set the material type of some specified elements
%   elem: the vector of elements specified
%   type: 1 for neo-hookean, 2 for linear, and more to come...

% assert((mtype == 1 || mtype == 2) || mtype == 3);

if length(Y) == 1
    obj.a = a;
    obj.b = b;
    
    obj.M = sparse(size(obj.M,1),size(obj.M,2));
    obj.mu = Y ./ ( 2 * (1 + P) );
    obj.lambda = ( Y .* P ) ./ ( (1 + P) .* (1 - 2 * P) );
    obj.rho = Rho;
    obj.Y = Y;
    obj.P = P;
    obj.Rho = Rho;
    obj.material_type = mtype;
    for t = 1:obj.NT
        %     assert(obj.elemMaterialType(elem(t)) == 0); % check the element haven't been initialized to any material type
        
        for e = obj.elem(t,:)
            for mi = (e-1)*3+1:e*3
                obj.M(mi,mi) = obj.M(mi,mi)+obj.W(t)/4 * Rho;
            end
        end
    end
    
    % after setting the material, calculate the rest_stiffness
%     obj.SetCurrentState(zeros(obj.Dim*obj.N,1));
%     obj.rest_stiffness = obj.StiffnessMatrix;
    
else
    assert(length(Y) == obj.NT)
    
    obj.a = a;
    obj.b = b;
    
    obj.M = sparse(size(obj.M,1),size(obj.M,2));
    
    obj.mu = zeros(obj.NT,1);
    obj.lambda = zeros(obj.NT,1);
    obj.rho = zeros(obj.NT,1);
    obj.Y = zeros(obj.NT,1);
    obj.P = zeros(obj.NT,1);
    obj.Rho = zeros(obj.NT,1);
    obj.material_type = zeros(obj.NT,1);
    
    for t = 1:obj.NT
        obj.mu(t) = Y(t) ./ ( 2 * (1 + P(t)) );
        obj.lambda(t) = ( Y(t) .* P(t) ) ./ ( (1 + P(t)) .* (1 - 2 * P(t)) );
        obj.rho(t) = Rho(t);
        obj.Y(t) = Y(t);
        obj.P(t) = P(t);
        obj.Rho(t) = Rho(t);
        obj.material_type(t) = mtype(t);
        for e = obj.elem(t,:)
            for mi = (e-1)*3+1:e*3
                obj.M(mi,mi) = obj.M(mi,mi)+obj.W(t)/4 * Rho(t);
            end
        end
    end
end
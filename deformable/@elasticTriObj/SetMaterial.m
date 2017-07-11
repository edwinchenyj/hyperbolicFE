function SetMaterial(obj, Y, P, Rho, elem, type)
% Set the material type of some specified elements
%   elem: the vector of elements specified
%   type: 1 for neo-hookean, 2 for linear, and more to come...
assert(max(elem) <= obj.NT);
assert(min(elem) > 0);
assert(type == 1 || type == 2);


obj.materialBlockCount = obj.materialBlockCount+1;
count = obj.materialBlockCount;



obj.materialBlock{count} = elem;
obj.materialType = [obj.materialType type];
obj.M = sparse(size(obj.M,1),size(obj.M,2));
for t = 1:length(elem)
%     assert(obj.elemMaterialType(elem(t)) == 0); % check the element haven't been initialized to any material type
    obj.elemMaterialType(elem(t)) = type;
    obj.mu(elem(t)) = Y ./ ( 2 * (1 + P) );
    obj.lambda(elem(t)) = ( Y .* P ) ./ ( (1 + P) .* (1 - 2 * P) );
    obj.rho(elem(t)) = Rho;
    % simple mass lumping
    % probably not right...
    for e = obj.elem(t,:)
        for mi = (e-1)*2+1:e*2
            obj.M(mi,mi) = obj.M(mi,mi)+obj.W(t)/3 * Rho;
        end
    end
end


end
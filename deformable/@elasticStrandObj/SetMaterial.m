function SetMaterial(obj, k, Rho, elem, type)
% Set the material type of some specified elements
%   elem: the vector of elements specified
%   type: 1 for neo-hookean, 2 for linear, and more to come...
%   % currently only linear works
assert(max(elem) <= obj.NS);
assert(min(elem) > 0);
assert(type == 1 | type == 2);


obj.materialBlockCount = obj.materialBlockCount+1;
count = obj.materialBlockCount;



obj.materialBlock{count} = elem;
obj.materialType = [obj.materialType type];

% obj.M = sparse(size(obj.M,1),size(obj.M,2));
for t = 1:length(elem)
    assert(obj.elemMaterialType(elem(t)) == 0); % check the element haven't been initialized to any material type
    obj.elemMaterialType(elem(t)) = type;
    obj.ks(elem(t)) = k;
    obj.rho(elem(t)) = Rho;
    for e = obj.elem(elem(t),:)
        for mi = (e-1)*3+1:e*3
            obj.M(mi,mi) = obj.M(mi,mi)+obj.W(t)/2 * Rho;
        end
    end
end


end
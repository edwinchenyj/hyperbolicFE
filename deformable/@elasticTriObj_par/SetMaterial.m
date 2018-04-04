function SetMaterial(obj, Y, P, Rho, type, a, b)
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
% simple mass lumping
% probably not right...
%     for e = obj.elem(t,:)
%         for mi = (e-1)*2+1:e*2
%             obj.M(mi,mi) = obj.M(mi,mi)+obj.W(t)/3 * Rho;
%         end
%     end
for t = 1:obj.NT
    % consistent mass from eq(31.27) in
    % http://kis.tu.kielce.pl/mo/COLORADO_FEM/colorado/IFEM.Ch31.pdf
    local_elem = obj.elem(t,:);
    i_elems = [ 2*local_elem-1; 2*local_elem];
    i_elems = i_elems(:);
    obj.M(i_elems,i_elems) = obj.M(i_elems,i_elems) +...
        Rho*obj.W(t)/12 *...
        [2 0 1 0 1 0;
        0 2 0 1 0 1;
        1 0 2 0 1 0;
        0 1 0 2 0 1;
        1 0 1 0 2 0;
        0 1 0 1 0 2];
    %     end
end


end
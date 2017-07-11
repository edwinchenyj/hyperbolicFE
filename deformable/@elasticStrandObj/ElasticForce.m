function f = ElasticForce(obj)
% compute the elastic force under the current deformation (3N by 1)
% can be calculated with current state as the sole input...

assert(obj.finalized);

    f = zeros(3*obj.N,1);
    
for s = 1:obj.NS

    ks = obj.ks(s);
    F = obj.F(3*(s-1)+1:3*s,:);
    l0 = obj.W(s);
    i = obj.elem(s, 1); j = obj.elem(s, 2);
    force = 1/2 * ks * (F'*F - 1) * F * l0;
    f(3*(i-1)+1:3*i) = f(3*(i-1)+1:3*i) + force;
    f(3*(j-1)+1:3*j) = f(3*(j-1)+1:3*j) - force;
end
end
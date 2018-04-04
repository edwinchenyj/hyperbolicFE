function df = ElasticForceDifferential(obj, dx)
% the matrix-free calculation of df = -K*dx, where K is the
% tangential stiffness matrix at the current deformation and dx
% is the small differential of x, which can be pictured as a
% small vector in the displacement field. This method can also
% be used to calculate the matrix-free product of K*v by
% substituting (dx) by (-v)

assert(obj.finalized)

% reshape the input dx into a nodal displacement field
dxNode = reshape(dx,3,obj.N)';

    df = zeros(3*obj.N,1);
for s = 1:obj.NS
    dx = dxNode(obj.elem(s,2),:)' - dxNode(obj.elem(s,1),:)';
    ks = obj.ks(s);
    F = obj.F(3*(s-1)+1:3*s,:);
    l0 = obj.W(s);
    dforce = 1/2 * ks * ((2*F' * dx) * F + (F' * F - 1)*dx);
    
    i = obj.elem(s,1); j = obj.elem(s,2);

    df(3*(i-1)+1:3*i) = df(3*(i-1)+1:3*i)+dforce(:,1);
    df(3*(j-1)+1:3*j) = df(3*(j-1)+1:3*j)-dforce(:,1);
end
end
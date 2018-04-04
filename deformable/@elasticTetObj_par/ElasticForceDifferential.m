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
for t = 1:obj.NT
    dDs = dxNode(obj.elem(t,:),:)' * obj.G;
    dF = dDs * obj.DmINV(3*(t-1)+1:3*t,:);
%     disp('dF= [')
%     disp(dF)    
%     disp('];')
    dP = obj.StressDifferential(t,dF);
%     disp('dP= [')
%     disp(dP)
%     disp('];')
    dH = -obj.W(t) *dP * obj.DmINV(3*(t-1)+1:3*t,:)';
    i = obj.elem(t,1); j = obj.elem(t,2); k = obj.elem(t,3); l = obj.elem(t,4);
    df(3*(i-1)+1:3*i) = df(3*(i-1)+1:3*i)+dH(:,1);
    df(3*(j-1)+1:3*j) = df(3*(j-1)+1:3*j)+dH(:,2);
    df(3*(k-1)+1:3*k) = df(3*(k-1)+1:3*k)+dH(:,3);
    df(3*(l-1)+1:3*l) = df(3*(l-1)+1:3*l) - dH(:,1) - dH(:,2) - dH(:,3);
end
end
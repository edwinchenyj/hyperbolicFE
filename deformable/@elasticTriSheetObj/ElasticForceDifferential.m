function df = ElasticForceDifferential(obj, dx)
% the matrix-free calculation of df = -K*dx, where K is the
% tangential stiffness matrix at the current deformation and dx
% is the small differential of x, which can be pictured as a
% small vector in the displacement field. This method can also
% be used to calculate the matrix-free product of K*v by
% substituting (dx) by (-v)

% reshape the input dx into a nodal displacement field
dxNode = reshape(dx,2,obj.N)';

    df = zeros(2*obj.N,1);
for t = 1:obj.NT
    dDs = dxNode(obj.elem(t,:),:)' * obj.G;
    dF = dDs * obj.DmINV(2*(t-1)+1:2*t,:);
%     disp('dF= [')
%     disp(dF)    
%     disp('];')
    dP = obj.StressDifferential(t,dF);
%     disp('dP= [')
%     disp(dP)
%     disp('];')
    dH = -obj.W(t) *dP * obj.DmINV(2*(t-1)+1:2*t,:)';
    i = obj.elem(t,1); j = obj.elem(t,2); k = obj.elem(t,3);
    df(2*(i-1)+1:2*i) = df(2*(i-1)+1:2*i)+dH(:,1);
    df(2*(j-1)+1:2*j) = df(2*(j-1)+1:2*j)+dH(:,2);
    df(2*(k-1)+1:2*k) = df(2*(k-1)+1:2*k) - dH(:,1) - dH(:,2);

end
end
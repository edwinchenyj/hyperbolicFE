function SetCurrentState(obj, Dx, v)
% set the deformation using Dx, the displacement field on the nodes, and v, 
% the velocity of the nodes. Notice Dx is not nodal position
Dx = Dx(:); v = v(:);
obj.x = obj.X + Dx;
obj.v = v;
obj.node = reshape(obj.x,3,obj.N)'; % update nodal positions
for i = 1:obj.NS
    S_node = obj.node(obj.elem(i,:),:); % element nodal position in the world space
    obj.Ds(3*(i-1)+1:3*i,:) = S_node' * obj.G;
    obj.F(3*(i-1)+1:3*i,:) = S_node' * obj.G / obj.W(i);
end
end
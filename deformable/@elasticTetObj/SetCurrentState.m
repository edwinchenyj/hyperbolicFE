function SetCurrentState(obj, Dx)
% set the deformation using Dx, the displacement field on the nodes, and v, 
% the velocity of the nodes. Notice Dx is not nodal position
Dx = Dx(:);
obj.x = obj.X + Dx;

obj.node = reshape(obj.x,3,obj.N)'; % update nodal positions
for i = 1:obj.NT
    T_node = obj.node(obj.elem(i,:),:); % element nodal position in the world space
    obj.Ds(3*(i-1)+1:3*i,:) = T_node' * obj.G;
    obj.F(3*(i-1)+1:3*i,:) = T_node' * obj.G * obj.DmINV(3*(i-1)+1:3*i,:);
    if obj.material_type == 1
        % only need FINV for neo-hookean
        obj.FINV(3*(i-1)+1:3*i,:) = inv(obj.F(3*(i-1)+1:3*i,:));
    end
end
end
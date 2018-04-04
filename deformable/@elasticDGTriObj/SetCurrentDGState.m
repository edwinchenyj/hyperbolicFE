function SetCurrentDGState(obj, DGDx)
% set the deformation using Dx, the displacement field on the nodes, and v, 
% the velocity of the nodes. Notice Dx is not nodal position
DGDx = DGDx(:);
obj.DGx = obj.DGX + DGDx;
obj.DGnode = reshape(obj.DGx,2,obj.DGN)'; % update nodal positions
for i = 1:obj.NT
    T_node = obj.DGnode(obj.DGelem(i,:),:); % element nodal position in the world space
    obj.Ds(2*(i-1)+1:2*i,:) = T_node' * obj.G;
    obj.F(2*(i-1)+1:2*i,:) = T_node' * obj.G * obj.DmINV(2*(i-1)+1:2*i,:);
    obj.FINV(2*(i-1)+1:2*i,:) = inv(obj.F(2*(i-1)+1:2*i,:));
end
obj.DGCalculateConst;

end
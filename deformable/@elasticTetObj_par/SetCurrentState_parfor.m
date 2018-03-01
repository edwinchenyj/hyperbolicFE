function SetCurrentState_parfor(obj, Dx)
% set the deformation using Dx, the displacement field on the nodes, and v, 
% the velocity of the nodes. Notice Dx is not nodal position
Dx = Dx(:);
obj.x = obj.X + Dx;

obj.node = reshape(obj.x,3,obj.N)'; % update nodal positions

node = obj.node;
elem = obj.elem;
G = obj.G;
Ds_flat = zeros(obj.NT,9);
F_flat = zeros(obj.NT,9);
FINV_flat = zeros(obj.NT,9);
DmINV_flat = reshape(reshape(obj.DmINV',1,[]),[],2)';
parfor i = 1:obj.NT
    
    T_node = node(elem(i,:),:); % element nodal position in the world space
    Ds_flat(i,:) = reshape(T_node' * G,1,[]);
    DmINV = reshape(DmINV_flat(i,:),3,[])';
    F = T_node' * G * DmINV;
    F_flat(i,:) = reshape(T_node' * G * DmINV,1,[]);
    FINV_flat(i,:) = inv(F);
end

end
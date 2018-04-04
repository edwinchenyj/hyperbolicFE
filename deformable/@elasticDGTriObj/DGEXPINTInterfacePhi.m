function phi = DGEXPINTInterfacePhi(obj,dt)

if isempty(obj.InterfacePhi)
    K_interface = obj.DGInterfaceStiffnessMatrix;
    A = zeros(size(K_interface,1)*2);
    A(end/2+1:end,1:end/2) = -dt*K_interface;
    [v,d] = eig(A);
    for id = 1:size(d,1)
        if abs(d(id,id)) < 1e-6
            d(id,id) = 1;
        else
            d(id,id) = (exp(d(id,id))-1)/d(id,id);
        end
    end
    obj.InterfacePhi = v*d*v';
    phi = obj.InterfacePhi;
else
    phi = obj.InterfacePhi;
end
end
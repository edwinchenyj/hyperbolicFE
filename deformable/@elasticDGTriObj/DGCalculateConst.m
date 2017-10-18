function DGCalculateConst(obj)
% Calculate b in
%       x = FX + b
% - this is constant is not used in linear CG elements, but it is used in
% linear DG elements
% - this is not constant in time, but constant in each element at each time
% step

obj.DGb = zeros(2*obj.NT,1);

for t = 1:obj.NT
    node = obj.DGnode(obj.DGelem(t,1),:)';
    nodeM = obj.DGnodeM(obj.DGelem(t,1),:)';
    F = obj.F(2*(t-1)+1:2*t,:);
    obj.DGb(2*(t-1)+1:2*t) = node - F*nodeM;
end


end
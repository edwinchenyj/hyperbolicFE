function mesh_quality(obj)
% calculate the ratio of area for the smallest and largest element
disp('mesh quality:')
AR = 0;
WR = 0;
Wmin = Inf;
Wmax = 0;
angle = Inf;
for i = 1:obj.NT
    T_node = obj.nodeM(obj.elem(i,:),:); % element nodal position in material space  (#nodes per elements by 2)
    n1 = T_node(1,:);
    n2 = T_node(2,:);
    n3 = T_node(3,:);
    e1 = norm(n1 - n2);
    e2 = norm(n2 - n3);
    e3 = norm(n3 - n1);
    E = [e1 e2 e3];
    AR = max([max(E)/min(E) AR]);
    Wmin = min([Wmin obj.W(i)]);
    Wmax = max([Wmax obj.W(i)]);
    angle1 = min([abs(acos((n1 - n2) * (n2 - n3)'/e1/e2)), pi - abs(acos((n1 - n2) * (n2 - n3)'/e1/e2))]);
    angle2 = min([abs(acos((n2 - n3) * (n3 - n1)'/e2/e3)), pi - abs(acos((n2 - n3) * (n3 - n1)'/e2/e3))]);
    angle3 = min([abs(acos((n3 - n1) * (n1 - n2)'/e3/e1)), pi - abs(acos((n3 - n1) * (n1 - n2)'/e3/e1))]);
    angle = min([angle angle1 angle2 angle3]);
end
disp('max aspect ratio')
disp(AR)
disp('area ratio')
disp(Wmax/Wmin)
disp('smallest angle')
disp(angle * 180/pi)
end
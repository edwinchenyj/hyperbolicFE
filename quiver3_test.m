cla
hold on
p1 = rand(3,1);
d1 = rand(3,2);

quiver3(p1(1),p1(2),p1(3),d1(1,1),d1(2,1),d1(3,1))
quiver3(p1(1),p1(2),p1(3),d1(1,2),d1(2,2),d1(3,2))

[q,r] = qr(d1);
r_new = [r [0 0 1]'];

d1_new = q*r_new
% cross(d1(:,1),d1(:,2))/norm(d1(:,1))/norm(d1(:,2));
d1_3 = d1_new(:,3);


quiver3(p1(1),p1(2),p1(3),d1_3(1),d1_3(2),d1_3(3))

d1_new'*d1_new


d2 = rand(3,2);

quiver3(p1(1),p1(2),p1(3),d2(1,1),d2(2,1),d2(3,1))
quiver3(p1(1),p1(2),p1(3),d2(1,2),d2(2,2),d2(3,2))

[q,r] = qr(d2);
r_new = [r [0 0 1]'];

d2_new = q*r_new
% cross(d2(:,1),d2(:,2))/norm(d2(:,1))/norm(d2(:,2));
d2_3 = d2_new(:,3);


quiver3(p1(1),p1(2),p1(3),d2_3(1),d2_3(2),d2_3(3))

d2_new'*d2_new

F = d1_new*inv(d2_new)
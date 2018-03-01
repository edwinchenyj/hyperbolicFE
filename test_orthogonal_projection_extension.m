n = 5;

K = rand(n);
K = 1/2 * (K + K');

M = full(gallery('tridiag',n,-1,2,-1));

[V,D] = eig(K,M);

V' * V
V * V'

V'*M*V
M*V*V'

R_s = diag([1.2,1.3]);
D_new1 = D;
D_new1(1:size(R_s,1),1:size(R_s,1)) = R_s * D(1:size(R_s,1),1:size(R_s,1));
K_new1 = M * V * D_new1 * (V') * M;

[V_1,D_1] = eig(K_new1,M)

V_s = V(:,1:size(R_s,1));
K_new2 = K + M * V_s * (-eye(size(R_s,1))) * (V_s') * K * V_s * (V_s') * M;

[V_2,D_2] = eig(K_new2,M)

K_new3 = M * V_s * (R_s) * (V_s') * K * V_s * (V_s') * M;

[V_3,D_3] = eig(K_new3,M)

K_new4 = K + M * V_s * (R_s-eye(size(R_s,1))) * (V_s') * K * V_s * (V_s') * M;


max(max(K_new1 - (K_new2 + K_new3)))

max(max(K_new1 - (K_new4)))


% 
% q = rand(n,1);
% 
% projected_q = V' * M * q
% back_extended_q = V * projected_q
% max(back_extended_q - q)
% 
% f = rand(n,1);
% 
% projected_f = V' * f;
% back_extended_f = M*V * projected_f;
% max(back_extended_f - f)
% 
% 
% % projected_K = V' * K * V;
% % back_extended_A = M*V*projected_K*V'*M;
% % max(max(back_extended_A - K))
% % 
% % v = rand(n,1);
% % 
% % projected_Mv = V' * M* v
% % back_extended_v = V * projected_Mv
% % max(back_extended_v - v)
% % 


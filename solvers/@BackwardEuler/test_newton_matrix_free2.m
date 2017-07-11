function test_newton_matrix_free2()
% 
% N = 10;
% 
% u0 = rand(N,1); % a random initial condition
% obj = BackwardEuler(u0');
% 
% J  = @(t, u) SPD(length(u)); % generate a SPD Jacobian
% JV = @(t, u, v) J(t, u) * v; % the matrix-vector product of the Jacobian with some vector (v)
% % notice the input also requires (t) & (u), which is required for the Jacobian although Jacobian may not be formed explicitly
% 
% obj.HJ = J;
% obj.HJV = JV;
% obj.KSMobj = CG();
% obj.KSMobj.SetHMV(@obj.ImhJ);
% 
% 
% 
% Jacobian = obj.HJ(t, u);
% h = 0.1;
% [D1, flag] = obj.KSMobj.solve();
% D2 = (obj.I - h * Jacobian) \ b;

end

function A = SPD(n)
% generating a SPD matrix;
    rng(1);
    A = rand(n);
    A = A + A';
    A = A + n*eye(n); % make sure it's diagonally dominating, for gershgorin circle theroem
end
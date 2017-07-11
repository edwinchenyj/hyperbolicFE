%%TEST3 testing wave equation with different sizes n = 10, 20 ,30, with
%% stiffness from 1 to 100000. print out error and
%%time at the end

m_max = [55];
p_max = [8];
i = [1:5]';
j = [1:1]';
k = 10.^i;
n = 10.*j;
% load('expm_exact.mat')

t0 = 0; tmax = 1;
tvals = tmax - t0;

funm_quad_time = zeros(length(i),1);
select_time = zeros(length(i),1);
backslash_time = zeros(length(i),1);
% choose parameters for the FUNM_QUAD restart algorithm
addpath('funm_quad')
param.function = 'exp';
param.restart_length = 70;          % each restart cycle consists of 70 Arnoldi iterations
param.max_restarts = 15;            % perform at most 15 restart cycles
param.tol = 1e-13;                  % tolerance for quadrature rule
param.hermitian = 0;                % the matrix A is Hermitian
param.V_full = 0;                   % set 1 if you need Krylov basis
param.H_full = 0;                   % do not store all Hessenberg matrices
param.exact = [];                   % exact solution. If not known set to []
param.stopping_accuracy = 1e-14;    % stopping accuracy
param.inner_product = @(a,b) b'*a;  % use standard Euclidean inner product
param.thick = [];                   % no implicit deflation is performed
param.min_decay = .95;              % we desire linear error reduction of rate < .95 
param.waitbar = 1;                  % show waitbar 
param.reorth_number = 0;            % #reorthogonalizations
param.truncation_length = inf;      % truncation length for Arnoldi 
param.verbose = 0;                  % print information about progress of algorithm

for n_index = 1:length(n)
    A = -gallery('poisson',n(n_index));
    b = linspace(-1,1,2 * n(n_index)^2)';
    A_tilde = sparse(2*size(A,1),2*size(A,1));
    A_tilde(1:(end/2), (end/2 + 1):end) = speye(size(A,1));
    
    for stiffness_index= 1:length(k)
        A_tilde((end/2 + 1):end, 1:(end/2)) = k(stiffness_index)*A;
        %         disp('selecting taylor degree')
        tic
        [X] = funm_quad(A_tilde,b,param);
        funm_quad_time(stiffness_index) = toc;
        
        % comparing the time for a single linear solve
        tic
        (tvals*A_tilde)\b;
        backslash_time(stiffness_index) = toc;
        
        %         fprintf('Relative differences between vectors from EXPM and EXPMV.\n')
        %         fprintf('with n = %d, s = %d, m = %d, k = %d.\n', [n, s(stiffness_index), m(stiffness_index), k(stiffness_index)]);
        %         fprintf('Should be of order %9.2e.\n', eps/2)
        Y = zeros(size(X));
        
        Y(:,1) = expm(full(tvals(1)*A_tilde))*b;
        err(stiffness_index) = norm(Y(:,1)-X(:,1),1)/norm(X(:,1),1);
        %         fprintf('%2.0f:  %9.2e\n', 1, err(stiffness_index));
        
    end
    fprintf('n = %d\n', n(n_index));
    fprintf('------------------\n')
    fprintf('stiffness k = %8.0f, err = %8.8e  funm_quad time = %8.8e, backslash time = %8.8e.\n', [k, err, funm_quad_time, backslash_time]')
    fprintf('funm_quad time / backslash time = %8.4e.\n',[funm_quad_time./backslash_time]);
end

%%TEST3 testing wave equation with different sizes n = 10, 20 ,30, with
%% stiffness from 1 to 100000. print out error and
%%time at the end

% load('expm_exact.mat')

j = [1:length(expm_exact)]';
n = 10.*j;
i = [1:size(expm_exact{1},2)]';
% i = 5;
k = 10.^i;

t0 = 0; tmax = 1;
tvals = tmax - t0;

expv_time = zeros(length(i),1);
select_time = zeros(length(i),1);
backslash_time = zeros(length(i),1);
err = zeros(length(i),1);
s = zeros(length(i),1);
m = zeros(length(i),1);
for n_index = 5:length(n)
    A = -gallery('poisson',n(n_index));
    b = linspace(-1,1,2 * n(n_index)^2)';
    A_tilde = sparse(2*size(A,1),2*size(A,1));
    A_tilde(1:(end/2), (end/2 + 1):end) = speye(size(A,1));
    
    for stiffness_index= 5:length(k)
        A_tilde((end/2 + 1):end, 1:(end/2)) = k(stiffness_index)*A;
        
        tic
        [X] = expv(tvals, A_tilde, b);
        expv_time(stiffness_index) = toc;
        
        % comparing the time for a single linear solve
        tic
        (tvals*A_tilde)\b;
        backslash_time(stiffness_index) = toc;
        
        %         fprintf('Relative differences between vectors from EXPM and EXPMV.\n')
        %         fprintf('with n = %d, s = %d, m = %d, k = %d.\n', [n, s(stiffness_index), m(stiffness_index), k(stiffness_index)]);
        %         fprintf('Should be of order %9.2e.\n', eps/2)
        Y = zeros(size(X));
        
        Y(:,1) = expm_exact{n_index}(:,stiffness_index);
        err(stiffness_index) = norm(Y(:,1)-X(:,1),1)/norm(X(:,1),1);
        %         fprintf('%2.0f:  %9.2e\n', 1, err(stiffness_index));
        
    end
    fprintf('n = %d\n', n(n_index));
    fprintf('------------------\n')
    fprintf('stiffness k = %8.0f, err = %8.8e, expv time = %8.8e, backslash time = %8.8e.\n', [k, err, expv_time, backslash_time]')
    fprintf('expv time / backslash time = %8.4e.\n',[expv_time./backslash_time]);
end

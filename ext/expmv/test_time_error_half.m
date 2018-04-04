%%TEST3 testing wave equation with different sizes n = 10, 20 ,30, with
%% stiffness from 1 to 100000. print out error and
%%time at the end

m_max = [55];
p_max = [8];
i = [1:5]';
k = 10.^i;

t0 = 0; tmax = 1;
tvals = tmax - t0;

expmv_time = zeros(length(i),1);
select_time = zeros(length(i),1);
backslash_time = zeros(length(i),1);
err = zeros(length(i),1);
s = zeros(length(i),1);
m = zeros(length(i),1);
%% size n = 10, small size, can still check error with expm from matlab
n = 10;
A = -gallery('poisson',n);
b = linspace(-1,1,2 * n^2)';
A_tilde = sparse(2*size(A,1),2*size(A,1));
A_tilde(1:(end/2), (end/2 + 1):end) = speye(size(A,1));

for stiffness_index= 1:length(k)
        A_tilde((end/2 + 1):end, 1:(end/2)) = k(stiffness_index)*A;
%         disp('selecting taylor degree')
        tic
        M = select_taylor_degree(A_tilde,b,m_max,p_max);
        
        select_time(stiffness_index) = toc;
        [X,s(stiffness_index), m(stiffness_index)] = expmv(tvals, A_tilde, b, M,'half');
        expmv_time(stiffness_index) = toc;
        
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
fprintf('n = %d\n', n);
fprintf('------------------\n')
fprintf('stiffness k = %8.0f, s = %2.0f, m = %2.0f, err = %8.8e, select time = %8.8e,  expmv time = %8.8e, backslash time = %8.8e.\n', [k, s, m, err, select_time, expmv_time, backslash_time]')


%% size n = 20, small size, can still check error with expm from matlab
n = 20;
A = -gallery('poisson',n);
b = linspace(-1,1,2 * n^2)';
A_tilde = sparse(2*size(A,1),2*size(A,1));
A_tilde(1:(end/2), (end/2 + 1):end) = speye(size(A,1));

for stiffness_index= 1:length(k)
        A_tilde((end/2 + 1):end, 1:(end/2)) = k(stiffness_index)*A;
%         disp('selecting taylor degree')
        tic
        M = select_taylor_degree(A_tilde,b,m_max,p_max);
        
        select_time(stiffness_index) = toc;
        [X,s(stiffness_index), m(stiffness_index)] = expmv(tvals, A_tilde, b, M,'half');
        expmv_time(stiffness_index) = toc;
        
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

fprintf('n = %d\n', n);
fprintf('------------------\n')
fprintf('stiffness k = %8.0f, s = %2.0f, m = %2.0f, err = %8.8e, select time = %8.8e,  expmv time = %8.8e, backslash time = %8.8e.\n', [k, s, m, err, select_time, expmv_time, backslash_time]')


%% size n = 30, small size, can still check error with expm from matlab
n = 30;
A = -gallery('poisson',n);
b = linspace(-1,1,2 * n^2)';
A_tilde = sparse(2*size(A,1),2*size(A,1));
A_tilde(1:(end/2), (end/2 + 1):end) = speye(size(A,1));

for stiffness_index= 1:length(k)
        A_tilde((end/2 + 1):end, 1:(end/2)) = k(stiffness_index)*A;
%         disp('selecting taylor degree')
        tic
        M = select_taylor_degree(A_tilde,b,m_max,p_max);
        
        select_time(stiffness_index) = toc;
        [X,s(stiffness_index), m(stiffness_index)] = expmv(tvals, A_tilde, b, M,'half');
        expmv_time(stiffness_index) = toc;
        
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

fprintf('n = %d\n', n);
fprintf('------------------\n')
fprintf('stiffness k = %8.0f, s = %2.0f, m = %2.0f, err = %8.8e, select time = %8.8e,  expmv time = %8.8e, backslash time = %8.8e.\n', [k, s, m, err, select_time, expmv_time, backslash_time]')

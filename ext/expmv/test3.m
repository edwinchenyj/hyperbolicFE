%TEST3 testing poisson with different size and stiffness, with m = 15, 30, 50.

%% size n = 10, small size, can still check error with expm from matlab
n = 10;
A = -gallery('poisson',n);
b = linspace(-1,1,2 * n^2)';
A_tilde = sparse(2*size(A,1),2*size(A,1));
A_tilde((end/2 + 1):end, 1:(end/2)) = A;
A_tilde(1:(end/2), (end/2 + 1):end) = speye(size(A,1));
t0 = 0; tmax = 1;
t = tmax - t0;
m_max = 55; 
p_max = 8;
disp('selecting taylor degree')
tic
M = select_taylor_degree(A_tilde,b,m_max,p_max);
toc
tvals = t;
tic
[X,s] = expmv(t, A_tilde, b, M);
toc

% comparing the time for a single linear solve
tic
(t*A_tilde)\b;
toc

fprintf('Relative differences between vectors from EXPM and EXPMV.\n')
fprintf('Should be of order %9.2e.\n', eps/2)
Y = zeros(size(X));
for i = 1:length(tvals)
    Y(:,i) = expm(full(tvals(i)*A_tilde))*b;
    fprintf('%2.0f:  %9.2e\n', i, norm(Y(:,i)-X(:,i),1)/norm(X(:,i),1) )
end

k = 10000;
A_tilde((end/2 + 1):end, 1:(end/2)) = k*A;
M = select_taylor_degree(A_tilde,b,m_max,p_max);
tic
[X,s] = expmv(t, A_tilde, b, M);
toc

% comparing the time for a single linear solve
tic
(t*A_tilde)\b;
toc

fprintf('Relative differences between vectors from EXPM and EXPMV.\n')
fprintf('Should be of order %9.2e.\n', eps/2)
Y = zeros(size(X));
for i = 1:length(tvals)
    Y(:,i) = expm(full(tvals(i)*A_tilde))*b;
    fprintf('%2.0f:  %9.2e\n', i, norm(Y(:,i)-X(:,i),1)/norm(X(:,i),1) )
end

%% size n = 20, small size, can still check error with expm from matlab
n = 20;
A = -gallery('poisson',n);
b = linspace(-1,1,2 * n^2)';
A_tilde = sparse(2*size(A,1),2*size(A,1));
A_tilde((end/2 + 1):end, 1:(end/2)) = A;
A_tilde(1:(end/2), (end/2 + 1):end) = speye(size(A,1));
t0 = 0; tmax = 1;
t = tmax - t0;
m_max = 55;
p_max = 8;
M = select_taylor_degree(A_tilde,b,m_max,p_max);
tic
[X] = expmv(t, A_tilde, b, M);
toc

% comparing the time for a single linear solve
tic
(t*A_tilde)\b;
time = toc

fprintf('Relative differences between vectors from EXPM and EXPMV.\n')
fprintf('Should be of order %9.2e.\n', eps/2)
Y = zeros(size(X));
for i = 1:length(tvals)
    Y(:,i) = expm(full(tvals(i)*A_tilde))*b;
    fprintf('%2.0f:  %9.2e\n', i, norm(Y(:,i)-X(:,i),1)/norm(X(:,i),1) )
end
k = 10000;
A_tilde((end/2 + 1):end, 1:(end/2)) = k*A;

M = select_taylor_degree(A_tilde,b,m_max,p_max);
tic
[X] = expmv(t, A_tilde, b, M);
toc

% comparing the time for a single linear solve
tic
(t*A_tilde)\b;
toc

fprintf('Relative differences between vectors from EXPM and EXPMV.\n')
fprintf('Should be of order %9.2e.\n', eps/2)
Y = zeros(size(X));
for i = 1:length(tvals)
    Y(:,i) = expm(full(tvals(i)*A_tilde))*b;
    fprintf('%2.0f:  %9.2e\n', i, norm(Y(:,i)-X(:,i),1)/norm(X(:,i),1) )
end

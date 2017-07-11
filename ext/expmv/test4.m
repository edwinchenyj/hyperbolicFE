%TEST   Simple test of EXPMV_TSPAN.
clear all
close all
n = 20;
A_bar = -1000*gallery('poisson',n);
b = linspace(-1,1,2 * n^2)';
A = sparse(2*size(A_bar,1),2*size(A_bar,1));
A((end/2 + 1):end, 1:(end/2)) = A_bar;
A(1:(end/2), (end/2 + 1):end) = speye(size(A_bar,1));
t0 = 0; tmax = 1;
q = 9;
t = tmax - t0;
[X] = expmv(t, A,b);
tvals = 1;

fprintf('Relative differences between vectors from EXPM and EXPMV_TSPAN.\n')
fprintf('Should be of order %9.2e.\n', eps/2)
Y = zeros(size(X));
for i = 1:length(tvals)
    Y(:,i) = expm(full(tvals(i)*A))*b;
    fprintf('%2.0f:  %9.2e\n', i, norm(Y(:,i)-X(:,i),1)/norm(X(:,i),1) )
end

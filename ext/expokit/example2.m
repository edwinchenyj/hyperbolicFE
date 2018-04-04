% MATLAB script illustrating the use of MEXPV.

disp('Loading the matrix ...');
n = 10;
A_bar = -gallery('poisson',n);
b = linspace(-1,1,2 * n^2)';
A = sparse(2*size(A_bar,1),2*size(A_bar,1));
A((end/2 + 1):end, 1:(end/2)) = A_bar;
A(1:(end/2), (end/2 + 1):end) = speye(size(A_bar,1));
t0 = 0; tmax = 1;
t = tmax - t0;

disp('Computing w = exp(A)e_1 ...');
[n,n] = size(A);
v = eye(n,1);
tic
[w,err] = expv(1,A*t,v);
toc
disp('w(1:10) =');
disp(w(1:10));

disp('err =');
disp(err)


tic
(A*t)\v;
toc
clear all
close all
n = 200;
c = 3;
d = 2;
e = 3;
A = gallery('tridiag',n,c,d,e);

[V,D] = eigs((A));

[D,ind] = sort(diag(D));
D  = diag(D);
D = sparse(D);
V = V(:,ind);

disp('original density of the sparse matrix')
nnz(A)/prod(size(A))

b = rand(size(A,1));

tic
for i = 1:1000
A*b;
end
disp('multiply the original matrix')
toc

% fill-in!

% modify the eigenvalues
D(1,1) = 0.1 * D(1,1);
D(2,2) = 0.1 * D(2,2);
D(3,3) = 0.1 * D(3,3);
A_processed = sparse((V*D*(V')));
A_processed(abs(A_processed) < 1e-8) = 0; 
% nnz(((V')*D*V))/prod([size((V')*D*V)])
disp('density of the modified sparse matrix')
nnz(A_processed)/prod([size(V*D*(V'))])

tic
for i = 1:1000
A_processed*b;
end
disp('multiply the modified matrix')
toc

tic
for i = 1:1000
V*D*(V')*b;
end
disp('multiply the modified matrix sequentially')
toc

clear all
close all
n = 20;
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
% A_processed = sparse(((V)*D*V'));
% A_processed(abs(A_processed) < 1e-8) = 0; 
% % nnz(((V')*D*V))/prod([size((V')*D*V)])
% nnz(A_processed)/prod([size(V*D*(V'))])
% 
% norm(A,1)
% % norm(A,2)
% norm(A,'inf')
% 
% norm(A_processed,1)
% % norm(A,2)
% norm(A_processed,'inf')
% 

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


% nnz(((V')*D*V))/prod([size((V')*D*V)])
% 
% norm(A,1)
% % norm(A,2)
% norm(A,'inf')
% 
% norm(A_processed,1)
% % norm(A,2)
% norm(A_processed,'inf')
% 
% 
% norm(A,1)
% normest(A,2)
% norm(A,'inf')

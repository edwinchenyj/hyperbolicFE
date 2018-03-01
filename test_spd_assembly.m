for it = 1:10000
Dim = 3;
A1 = rand(Dim);
A1 = A1 + A1';
[V1,D1] = eig(A1);

A2 = rand(Dim);
A2 = A2 + A2';
[V2,D2] = eig(A2);

for i = 1:Dim
    if D1(i,i)<0
        D1(i,i) = 1e-6;
    end
    
    if D2(i,i)<0
        D2(i,i) = 1e-6;
    end
end

A1 = V1 * D1 * inv(V1);
A2 = V2 * D2 * inv(V2);
A = zeros(Dim+1);
A(1:end-1,1:end-1) = A1;
A(2:end,2:end) = A(2:end,2:end) + A2;
disp('eig')
disp(eig(A))
% any(eig(A) < 0)
assert(any(eig(A)>0))
end
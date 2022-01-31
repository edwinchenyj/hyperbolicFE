
msinc(zeros(5))

function out = msinc(A)
n = size(A,1);

[U,D] = eig(A);

for j = 1:n
    if norm(D(j,j)) > 1e-8
        D(j,j) = sin(D(j,j))/D(j,j);
    else
        D(j,j) = 1;
    end
end

out = real(U*D/(U));
end

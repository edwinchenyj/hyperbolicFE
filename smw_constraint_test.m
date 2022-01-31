n = 10;

constraint_index = 5;

A = diag(1:n);
b = rand(n,1);

Y = rand(n,1);
Z = rand(n,1);

A(constraint_index,constraint_index) = 1;
Y(constraint_index) = 0;
Z(constraint_index) = 0;


A_new = A + Y * (Z');


result1 = A_new\b;

result2 = A\b - (A\Y)*((speye(1) + (Z')*(A\Y))\(Z')*(A\b))

result1 - result2

b(constraint_index) - result1(constraint_index)
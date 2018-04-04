n = 50;
A = gallery('poisson',n);

b = rand(n^2,1);
tic
B = inv(A);
B*b;
toc

tic
A\b;
A\b;
A\b;
A\b;
A\b;
toc
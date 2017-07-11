function test_rosenbrock2()

close all

% testing with a simple ode y' = -y;
f = @(t, y) (-y);
Jv = @(t, u) -1;
g = @(t, u) 0;
y0 = 5;
obj = Rosenbrock2(f, y0, Jv, false);
h = .1;
t0 = 0;
obj.T = t0;
t = t0;
T = [t0];

tic
for i = 1:10
obj.step(h);
t = t + h;
T = [T; t];
end
toc

exact = @(x) y0*exp(-x);
plot(T, obj.state, 'x', T, exact(T))
legend('Rosenbrock','exact')

% testing with a simple system y' = A*y;
n = 100;
A = -rand(n)/n;
f1 = @(t, y) (A*y);
Jv = @(t, u) A;
y01 = rand(1,n);
obj2 = Rosenbrock2(f1, y01, Jv, false);
t0 = 0;
obj2.T = t0;
t = t0;
T = [t0];

tic
for i = 1:1000
obj2.step(h);
t = t + h;
T = [T; t];
end
toc

exact = @(x) expm(A*x)*y01';
exact_sol = [];
for i = T'
    exact_sol = [exact_sol exact(i)];
end
figure
plot(T, obj2.state, 'x', T, exact_sol')

error = norm(obj2.state(end,:) - exact_sol(:,end)')


%% using KSM to approximate matrix exponentials
% testing with a simple ode y' = -y;
f = @(t, y) (-y);
Jv = @(t, u) -1;
y0 = 5;
obj = Rosenbrock2(f, y0, Jv, true);
t0 = 0;
obj.T = t0;
t = t0;
T = [t0];

tic
for i = 1:10
obj.step(h);
t = t + h;
T = [T; t];
end
toc

figure
exact = @(x) y0*exp(-x);
plot(T, obj.state, 'x', T, exact(T))
legend('Rosenbrock','exact')

% testing with a simple system y' = A*y;
f1 = @(t, y) (A*y);
Jv = @(t, u) A;
obj2 = Rosenbrock2(f1, y01, Jv, true);
t0 = 0;
obj2.T = t0;
t = t0;
T = [t0];

tic
for i = 1:10
obj2.step(h);
t = t + h;
T = [T; t];
end
toc
exact = @(x) expm(A*x)*y01';
exact_sol = [];
for i = T'
    exact_sol = [exact_sol exact(i)];
end
figure
plot(T, obj2.state, 'x',T, exact_sol')

error = norm(obj2.state(end,:) - exact_sol(:,end)')

clear all

 
end


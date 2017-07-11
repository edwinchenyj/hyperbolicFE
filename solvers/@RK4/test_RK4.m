function test_RK4()
% testing with a simple ode y' = y;
f = @(t, y) y;
y0 = 5;
obj = RK4(f,y0);
h = 0.005;
t0 = 0;
obj.T = t0;
t = t0;
T = [t0];
for i = 1:1000
obj.step(h);
t = t + h;
T = [T; t];
end

exact = @(x) y0*exp(x);
plot(T, obj.state, T, exact(T))

% testing with a simple system y' = A*y;
A = rand(5);
f = @(t, y) (A*(y'))';
y0 = rand(1,5);
obj2 = RK4(f,y0);
h = 0.0005;
t0 = 0;
obj2.T = t0;
t = t0;
T = [t0];
for i = 1:1000
obj2.step(h);
t = t + h;
T = [T; t];
end

exact = @(x) expm(A*x)*y0';
exact_sol = [];
for i = T'
    exact_sol = [exact_sol exact(i)];
end
plot(T, obj2.state, T, exact_sol')


clear all;
end
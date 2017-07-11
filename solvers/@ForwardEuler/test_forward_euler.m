function test_forward_euler()
% testing with a simple ode y' = y;
f = @(t, y) y;
y0 = 5;
FE = ForwardEuler(f,y0);
h = 0.005;
t0 = 0;
FE.T = t0;
t = t0;
T = [t0];
for i = 1:1000
FE.step(h);
t = t + h;
T = [T; t];
end

exact = @(x) y0*exp(x);
plot(T, FE.state, T, exact(T))

% testing with a simple system y' = A*y;
A = rand(5);
f = @(t, y) (A*(y'))';
y0 = rand(1,5);
FE2 = ForwardEuler(f,y0);
h = 0.0005;
t0 = 0;
FE2.T = t0;
t = t0;
T = [t0];
for i = 1:1000
FE2.step(h);
t = t + h;
T = [T; t];
end

exact = @(x) expm(A*x)*y0';
exact_sol = [];
for i = T'
    exact_sol = [exact_sol exact(i)];
end
plot(T, FE2.state, T, exact_sol')

end
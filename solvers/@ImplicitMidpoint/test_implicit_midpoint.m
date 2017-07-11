function test_implicit_midpoint()
clear all
close all

% testing with a simple ode y' = -y;
f = @(t, y) (-y);
Jv = @(t, u) -1;
g = @(t, u) 0;
y0 = 5;
obj = ImplicitMidpoint(f, y0,false, Jv);
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

exact = @(x) y0*exp(-x);
plot(T, obj.state, T, exact(T))
legend('Midpoint','exact')

% testing with a simple system y' = A*y;

A = -rand(5);
f1 = @(t, y) (A*y);
Jv = @(t, u) A;
y01 = rand(1,5);
obj2 = ImplicitMidpoint(f1, y01, false, Jv);
h = 0.005;
t0 = 0;
obj2.T = t0;
t = t0;
T = [t0];
for i = 1:1000
    obj2.step(h);
    t = t + h;
    T = [T; t];
end

exact = @(x) expm(A*x)*y01';
exact_sol = [];
for i = T'
    exact_sol = [exact_sol exact(i)];
end
figure
plot(T, obj2.state, T, exact_sol')



%% matrix free
% testing with a simple ode y' = -y;
Jv = @(t, v) -v;
Jou = @(t, u) 1;
y0 = 5;
obj = ImplicitMidpoint(f,y0,true, Jv, Jou);
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

exact = @(x) y0*exp(-x);
figure
plot(T, obj.state, T, exact(T))

% testing with a simple system y' = A*y;
Jv = @(t, v) A*v;
Jou = @(t, u) 1;
obj2 = ImplicitMidpoint(f1, y01, true, Jv, Jou);
h = 0.005;
t0 = 0;
obj2.T = t0;
t = t0;
T = [t0];
for i = 1:1000
    obj2.step(h);
    t = t + h;
    T = [T; t];
end

exact = @(x) expm(A*x)*y01';
exact_sol = [];
for i = T'
    exact_sol = [exact_sol exact(i)];
end
figure
plot(T, obj2.state, T, exact_sol')


clear all


end


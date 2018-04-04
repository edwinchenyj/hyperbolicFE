function test_newton
f = @(x) x^2-2;
df = @(x) 2*x;
newton = Newton(f,1,df,false);
newton.solve
end
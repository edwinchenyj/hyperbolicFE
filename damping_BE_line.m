w = logspace(1,5,100);
d = 0.05*w;

t = 1e-1;
damping2 = 1/t * log(1+(t*w).^2);

plot(w,d,w,damping2)
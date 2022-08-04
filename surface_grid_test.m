close all

[X,Y] = meshgrid(1:0.5:100,1:200);
Z = sin(X) + cos(Y);
surface(X,Y,Z)
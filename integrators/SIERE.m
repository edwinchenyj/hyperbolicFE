function u_new = SIERE( dt, u, obj,reduced_n, varargin)
% inputs:
%   dt: step size
%    u: current state
%  obj: the elastic obj
%  varargin: optional input for constraints
% notice sum(nnz(~constraint_indices)) = length(u)/2

if nargin == 4
    constraint_indices = false(size(u,1)/2,1);
elseif nargin == 5
    constraint_indices = varargin{1};
end

indLogical = ~constraint_indices;

x = u(1:end/2);
v = u(end/2+1:end);
x = x(indLogical);
v = v(indLogical);

n = size(x,1);

K = obj.StiffnessMatrix;
M = obj.M;
% Eforce = obj.ElasticForce;

M = M(indLogical,indLogical);
K = K(indLogical,indLogical);
ElasticForce = obj.ElasticForce;
ElasticForce = ElasticForce(indLogical);

ExternalForce = M * obj.externalGravity(indLogical);

f = ElasticForce + ExternalForce;


[Us,D] = eigs(K,M,reduced_n,'smallestabs');

vG = Us*(Us')*M*v;
vH = v - vG;
fG = M*Us*(Us')*f;
fH = f - fG;

H = [vH; M\fH];

JGr = [zeros(reduced_n) eye(reduced_n); -D zeros(reduced_n)];
Gr = [Us'*M*v; Us'*f];

p = phi(dt * JGr);

delta = dt*H + dt * [Us zeros(n,reduced_n); zeros(n,reduced_n) Us] * p * Gr;

J = [zeros(n), speye(n); -M\K zeros(n)];
JG = [zeros(n) Us*(Us')*M; -Us*(Us')*K*Us*(Us')*M zeros(n)];
JH = J - JG;
x0 = (speye(2*n)-dt*JH)\delta;

x_new = u(1:end/2);
v_new = u(end/2+1:end);
x_new(indLogical) = x_new(indLogical) + x0(1:end/2);
v_new(indLogical) = v_new(indLogical) + x0(end/2+1:end);

obj.x = obj.X +  x_new;


obj.SetCurrentState(obj.x - obj.X);

u_new = [x_new;v_new];

end

function out = phi(A)
n = size(A,1);

[U,D] = eig(A);

for j = 1:n
    if norm(D(j,j)) > 1e-8
        temp = exp(D(j,j)) - 1;
        temp = temp/D(j,j);
        D(j,j) = temp;
    else
        D(j,j) = 1;
    end
end

out = real(U*D/(U));
end
function u_new = SIERE( dt, u, obj, reduced_n, constraint_indices, recompute, varargin)
% inputs:
%   dt: step size
%    u: current state
%  obj: the elastic obj
%  varargin: optional input for constraints
% notice sum(nnz(~constraint_indices)) = length(u)/2


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

if recompute
[Us,D] = eigs(K,M,reduced_n,'smallestabs');
obj.eig_modes = Us;
obj.eig_values = D;
else
Us = obj.eig_modes;
D = obj.eig_values;
end

obj.ritz_errors = vecnorm(M\K*Us - Us*D);

vG = Us*(Us')*M*v;
vH = v - vG;
fG = M*Us*(Us')*f;
fH = f - fG;

H = [vH; M\fH];

JGr = [zeros(reduced_n) eye(reduced_n); -D zeros(reduced_n)];
Gr = [Us'*M*v; Us'*f];

p = phi(dt * JGr);
%p_simple = phi_simple(dt,D);

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
% 
% function out = phi_simple(dt,D)
% evs = sqrt(diag(D)) * 1i * dt;
% evs = [evs; -evs];
% eVecs = zeros(size(D,1)*2);
% 
% for j = 1:length(evs)/2
%     eVecs(j,j) = 1;
%     eVecs(j+length(evs)/2,j) = evs(j);
%     eVecs(j,j+length(evs)/2) = -1;
%     eVecs(j+length(evs)/2,j+length(evs)/2) = evs(j);
% end
% evs = (exp(evs) -1)./evs;
% 
% 
% end
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

if ~isempty(obj.eig_values)
    obj.ritz_errors = vecnorm(M\K*obj.eig_modes - obj.eig_modes*obj.eig_values);
end

if norm(obj.ritz_errors) > 1e4
    disp('ritz error norm:')
    disp(norm(obj.ritz_errors))
    if sum(~indLogical) < 6
        [Us,D] = eigs(K,M,reduced_n + 6,'smallestabs');
        obj.eig_modes = Us(:,end-reduced_n+1:end);
        obj.eig_values = D(end-reduced_n+1:end,end-reduced_n+1:end);
        
    else
        [Us,D] = eigs(K,M,reduced_n,'smallestabs');
        obj.eig_modes = Us;
        obj.eig_values = D;
    end
end
Us = obj.eig_modes;
D = obj.eig_values;



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

disp('cond siere:')
disp(cond((speye(2*n)-dt*JH)))
disp('cond si:')
disp(cond((speye(2*n)-dt*J)))

disp('svd ratio siere:')
ev_large_siere = svds((speye(2*n)-dt*JH),1,'largest');
ev_small_siere = svds((speye(2*n)-dt*JH),1,'smallest');
disp(abs(ev_large_siere)/abs(ev_small_siere))
disp('svd ratio si:')
ev_large_si = svds((speye(2*n)-dt*J),1,'largest');
ev_small_si = svds((speye(2*n)-dt*J),1,'smallest');
disp(abs(ev_large_si)/abs(ev_small_si))


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
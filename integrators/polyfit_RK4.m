function u_new = polyfit_RK4( dt, u, obj, varargin)

if nargin == 3
    constraint_indices = false(size(u,1)/2,1);
elseif nargin == 4
    constraint_indices = varargin{1};
end


indLogical = ~constraint_indices;

Dx0 = u(1:end/2);
v0 = u(end/2+1 : end);
% k1
Dx = u(1:end/2);
v = u(end/2 + 1:end);

K = obj.StiffnessMatrix;
Mass = obj.M;
Mass = Mass(indLogical,indLogical);
K = K(indLogical,indLogical);

Minv_K_new = polyvalm([obj.polyfit_p],Mass\K);

Eforce = obj.ElasticForce;
Eforce = Mass*Minv_K_new*(K\Eforce(indLogical));
fExternal = Mass * obj.externalGravity(indLogical);


%         f1 = Eforce + fExternal + B*v(indLogical); % from column to row
f1 = Eforce + fExternal;

k1 = dt * [v(indLogical); Mass\f1];

% k2
Dx(indLogical) = Dx0(indLogical) + k1(1:end/2)/2;
v(indLogical) = v0(indLogical) + k1(end/2+1 : end)/2;

obj.SetCurrentState(Dx);
K = obj.StiffnessMatrix;
Mass = obj.M;
Mass = Mass(indLogical,indLogical);
K = K(indLogical,indLogical);

Minv_K_new = polyvalm([obj.polyfit_p],Mass\K);

Eforce = obj.ElasticForce;
Eforce = Mass*Minv_K_new*(K\Eforce(indLogical));
fExternal = Mass * obj.externalGravity(indLogical);

%         f2 = Eforce + fExternal + B*v(indLogical); % from column to row
f2 = Eforce + fExternal;

k2 = dt * [v(indLogical); Mass\f2];

% k3
Dx(indLogical) = Dx0(indLogical) + k2(1:end/2)/2;
v(indLogical) = v0(indLogical) + k2(end/2+1 : end)/2;

obj.SetCurrentState(Dx);
K = obj.StiffnessMatrix;
Mass = obj.M;
Mass = Mass(indLogical,indLogical);
K = K(indLogical,indLogical);

Minv_K_new = polyvalm([obj.polyfit_p],Mass\K);

Eforce = obj.ElasticForce;
Eforce = Mass*Minv_K_new*(K\Eforce(indLogical));
fExternal = Mass * obj.externalGravity(indLogical);

f3 = Eforce + fExternal;

k3 = dt * [v(indLogical); Mass\f3];

% k4
Dx(indLogical) = Dx0(indLogical) + k3(1:end/2);
v(indLogical) = v0(indLogical) + k3(end/2+1 : end);

obj.SetCurrentState(Dx);
K = obj.StiffnessMatrix;
Mass = obj.M;
Mass = Mass(indLogical,indLogical);
K = K(indLogical,indLogical);

Minv_K_new = polyvalm([obj.polyfit_p],Mass\K);

Eforce = obj.ElasticForce;
Eforce = Mass*Minv_K_new*(K\Eforce(indLogical));
fExternal = Mass * obj.externalGravity(indLogical);

f4 = Eforce + fExternal;

k4 = dt * [v(indLogical); Mass\f4];

% RK4
%         positions(indLogical) = positionsM(indLogical) + u(1:end/2) + k1(1:end/2)/6 + k2(1:end/2)/3 + k3(1:end/2)/3 + k4(1:end/2)/6;
%         v(indLogical) = u(end/2 + 1 : end) + k1(end/2 + 1 : end)/6 + k2(end/2 + 1 : end)/3 + k3(end/2 + 1 : end)/3 + k4(end/2 + 1 : end)/6;
Dx = u(1:end/2);
v = u(end/2+1 : end);
Dx(indLogical) = Dx(indLogical) +  k1(1:end/2)/6 + k2(1:end/2)/3 + k3(1:end/2)/3 + k4(1:end/2)/6;
v(indLogical) = v(indLogical) + k1(end/2 + 1 : end)/6 + k2(end/2 + 1 : end)/3 + k3(end/2 + 1 : end)/3 + k4(end/2 + 1 : end)/6;

u_new = [Dx; v];
end
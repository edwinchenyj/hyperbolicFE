function u_new = IMIMEX( dt, u, obj, varargin)
% inputs:
%   dt: step size
%    u: current state
%  obj: the elastic obj
%  varargin: optional input for constraints
% notice sum(nnz(~constraint_indices)) = length(u)/2

if nargin == 3
    constraint_indices = false(size(u,1)/2,1);
    MaxIT = 40;
elseif nargin == 4
    constraint_indices = varargin{1};
    MaxIT = 40;
elseif nargin == 5
    constraint_indices = varargin{1};
    MaxIT = varargin{2};
end

indLogical = ~constraint_indices;

it = 0;
dq = u(1:end/2);
v = u(end/2+1:end);
u_new = u;
v_new = u_new(end/2+1:end);

obj.x = obj.X + dq + 1/4 * dt * (v + v_new);

obj.SetCurrentState(obj.x - obj.X); % update deformation gradient
K_mid = obj.DGElementStiffnessMatrix;
Mass = obj.M;
Eforce_mid = obj.ElasticForce;
Ki = obj.DGInterfaceStiffnessMatrix;
Mass = Mass(indLogical,indLogical);
K_mid = K_mid(indLogical,indLogical);
Ki = Ki(indLogical,indLogical);
B = -obj.a * Mass - obj.b * Ki;

Eforce_mid = Eforce_mid(indLogical);

fExternal = Mass * obj.externalGravity(indLogical);

f_mid = Eforce_mid + fExternal + B*1/2*(v(indLogical) + v_new(indLogical));

N = obj.N;
nFixed = sum(constraint_indices)/obj.Dim;

residual0 = (dt * (Mass\f_mid))' * (dt * (Mass\f_mid));
Dv = -(speye(2*(N-nFixed)) + 1/4* dt*dt*(Mass\K_mid) - 1/2 * dt*(Mass\B))\(v_new(indLogical) - v(indLogical) - dt * (Mass\f_mid));
v_new(indLogical) = v_new(indLogical) + Dv;

residual = (v_new(indLogical) - v(indLogical) - dt * (Mass\f_mid))' * (v_new(indLogical) - v(indLogical) - dt * (Mass\f_mid));

it = it + 1;

u_new(1:end/2) = dq + 1/2 * dt * (v + v_new);
u_new(end/2+1:end) = v_new;

while (Dv'*Dv > 1e-12) && (residual > 1e-12)
    
    v_new = u_new(end/2+1:end);
    obj.x = obj.X + dq + 1/4*dt*(v + v_new);
    
    obj.SetCurrentState(obj.x - obj.X);
    K_mid = obj.DGElementStiffnessMatrix;
    Mass = obj.M;
    Eforce_mid = obj.ElasticForce;
    
    Mass = Mass(indLogical,indLogical);
    K_mid = K_mid(indLogical,indLogical);
    
    B = -obj.a * Mass - obj.b * Ki;
    
    Eforce_mid = Eforce_mid(indLogical);
    
    fExternal = Mass * obj.externalGravity(indLogical);
    
    f_mid = Eforce_mid + fExternal + B*1/2*(v(indLogical)+v_new(indLogical));
    
    Dv = -(speye(2*(N-nFixed)) + 1/4* dt*dt*(Mass\K_mid) - 1/2 * dt*(Mass\B))\(v_new(indLogical) - v(indLogical) - dt * (Mass\f_mid));
    v_new(indLogical) = v_new(indLogical) + Dv;
    
    residual = (v_new(indLogical) - v(indLogical) - dt * (Mass\f_mid))' * (v_new(indLogical) - v(indLogical) - dt * (Mass\f_mid));
    it = it + 1;
    
    u_new(1:end/2) = dq + 1/2 * dt * (v + v_new);
    u_new(end/2+1:end) = v_new;
    
    if (it > 3 && residual > residual0) || it == MaxIT
        disp('local substep required')
        if nargin > 3
            u_half = IMIMEX(dt/2, u, obj, varargin{:});
        else
            u_half = IMIMEX(dt/2, u, obj);
        end
        v = u_half(end/2 + 1:end);
        dq = u_half(1:end/2);
        
        obj.x = obj.X + dq;
        
        u_half = [obj.x-obj.X; v];
        
        if nargin > 3
            u_new = IMIMEX(dt/2, u_half, obj, varargin{:});
        else
            u_new = IMIMEX(dt/2, u_half, obj);
        end
        break;
    end
    
    if it == MaxIT
        error('Newton iteration not converging in IMIMEX')
    end
    
end

end
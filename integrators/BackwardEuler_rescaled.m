function u_new = BackwardEuler_rescaled( dt, u, obj, varargin)
% inputs:
%   dt: step size
%    u: current state
%  obj: the elastic obj
%  varargin: optional input for constraints and max iteration count
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
Dim = obj.Dim;

it = 0;
dq = u(1:end/2);
v = u(end/2+1:end);
u_new = u;
v_new = u_new(end/2+1:end);

obj.x = obj.X + dq + dt * (v_new);

obj.SetCurrentState(obj.x - obj.X); % update deformation gradient

Mass = obj.M;
Mass = Mass(indLogical,indLogical);
[K_end, Eforce_end] = obj.eigModification_nonlinear_subspace;


fExternal = Mass * obj.externalGravity(indLogical);

B = -obj.a * Mass - obj.b * K_end;

f_mid = Eforce_end + fExternal + B*v_new(indLogical);

N = obj.N;
nFixed = sum(constraint_indices)/obj.Dim;

residual0 = (dt * (Mass\f_mid))' * (dt * (Mass\f_mid));
Dv = -(speye(Dim*(N-nFixed)) + dt*dt*(Mass\K_end) - dt*(Mass\B))\(v_new(indLogical) - v(indLogical) - dt * (Mass\f_mid));
v_new(indLogical) = v_new(indLogical) + Dv;

residual = (v_new(indLogical) - v(indLogical) - dt * (Mass\f_mid))' * (v_new(indLogical) - v(indLogical) - dt * (Mass\f_mid));

it = it + 1;

u_new(1:end/2) = dq + dt * (v_new);
u_new(end/2+1:end) = v_new;

while (Dv'*Dv > 1e-12) && (residual > 1e-12)
    
    v_new = u_new(end/2+1:end);
    obj.x = obj.X + dq + dt*(v_new);
    
    obj.SetCurrentState(obj.x - obj.X);

    Mass = obj.M;
    Mass = Mass(indLogical,indLogical);
    [K_end, Eforce_end] = obj.eigModification_nonlinear_subspace;

    B = -obj.a * Mass - obj.b * K_end;
    
    fExternal = Mass * obj.externalGravity(indLogical);
    
    f_mid = Eforce_end + fExternal + B*(v_new(indLogical));
    
    Dv = -(speye(Dim*(N-nFixed)) + dt*dt*(Mass\K_end) - dt*(Mass\B))\(v_new(indLogical) - v(indLogical) - dt * (Mass\f_mid));
    v_new(indLogical) = v_new(indLogical) + Dv;
    
    residual = (v_new(indLogical) - v(indLogical) - dt * (Mass\f_mid))' * (v_new(indLogical) - v(indLogical) - dt * (Mass\f_mid));
    it = it + 1;
    
    u_new(1:end/2) = dq + dt * (v_new);
    u_new(end/2+1:end) = v_new;
    
    if (it > 3 && residual > residual0) || it == MaxIT
        disp('local substep required')
        if nargin > 3
            u_half = BackwardEuler_rescaled(dt/2, u, obj, varargin{1});
        else
            u_half = BackwardEuler_rescaled(dt/2, u, obj);
        end
        v = u_half(end/2 + 1:end);
        dq = u_half(1:end/2);
        
        obj.x = obj.X + dq;
        
        u_half = [obj.x-obj.X; v];
        
        if nargin > 3
            u_new = BackwardEuler_rescaled(dt/2, u_half, obj, varargin{1});
        else
            u_new = BackwardEuler_rescaled(dt/2, u_half, obj);
        end
        break;
    end
    
    if it == MaxIT
        error('Newton iteration not converging in BE')
    end
    
end

end
function  [u_new, residual, step_length] = ImplicitMid_rescaled_line_search( dt, u, obj, varargin)
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


dq = u(1:end/2);
v = u(end/2+1:end);
u_new = u;
v_new = u_new(end/2+1:end);
Mass = obj.M;
Mass = Mass(indLogical,indLogical);


Dv = Inf;
residual = Inf;

it_outer = 0;
while and(residual > 1e-4,  norm(Dv) > 1e-4)
    it_outer = it_outer + 1;
    if it_outer == 40
        error('quasi-newton more than 20 iterations.');
    end
    
    it = 0;
    obj.x = obj.X + dq + 1/4 * dt * (v + v_new);
    obj.SetCurrentState(obj.x - obj.X); % update deformation gradient
    
    [K_mid, Eforce_mid] = obj.eigModification_nonlinear_subspace;
    % else
    % assert(all(obj.x == obj.X + dq + 1/4 * dt * (v + v_new)));
    % K_mid = obj.K_prev;
    % Eforce_mid = obj.f_prev;
    % end
    
    B = -obj.a * Mass - obj.b * K_mid;
    
    fExternal = Mass * obj.externalGravity(indLogical);
    
    f_mid = Eforce_mid + fExternal + B*1/2*(v(indLogical) + v_new(indLogical));
    
    N = obj.N;
    nFixed = sum(constraint_indices)/obj.Dim;
    
    residual_old = (1/2) * (dt * (Mass\f_mid))' * (dt * (Mass\f_mid));
    
    % calculate descent direction
    Dv = -(speye(obj.Dim*(N-nFixed)) + 1/4* dt*dt*(Mass\K_mid) - 1/2 * dt*(Mass\B))\(v_new(indLogical) - v(indLogical) - dt * (Mass\f_mid));
    
    step_size = 1;
    v_new_temp = v_new;
    c1 = 1e-4;
    c2 = 0.9;
    
    v_new_temp(indLogical) = v_new(indLogical) + step_size*Dv;
    dq_new_temp = dq + (1/4) * dt * (v_new_temp + v);
    obj.SetCurrentState(dq_new_temp);
    
    Mass = obj.M;
    Mass = Mass(indLogical,indLogical);
    Eforce_mid = obj.eigModification_nonlinear_subspace_f_only;
    
    f_mid = Eforce_mid + fExternal + B*1/2*(v(indLogical) + v_new_temp(indLogical));
    
    % merit function
    residual = (1/2) * (v_new_temp(indLogical) - v(indLogical) - dt * (Mass\f_mid))' * (v_new_temp(indLogical) - v(indLogical) - dt * (Mass\f_mid));
    
    while (residual > residual_old + c1 * step_size * residual * 2)
        residual_old = residual;
        it = it + 1;
        if it == 40
            error('line search more than 40 iterations.');
        end
        step_size = 0.9 * step_size;
        v_new_temp(indLogical) = v_new(indLogical) + step_size*Dv;
        dq_new_temp = dq + (1/4) * dt * (v_new_temp + v);
        obj.SetCurrentState(dq_new_temp);
        
        Mass = obj.M;
        Mass = Mass(indLogical,indLogical);
        Eforce_mid = obj.eigModification_nonlinear_subspace_f_only;
        
        f_mid = Eforce_mid + fExternal + B*1/2*(v(indLogical) + v_new_temp(indLogical));
        
        % merit function
        residual = (1/2) * (v_new_temp(indLogical) - v(indLogical) - dt * (Mass\f_mid))' * (v_new_temp(indLogical) - v(indLogical) - dt * (Mass\f_mid));
        
    end
    
    v_new = v_new_temp;
end

u_new(1:end/2) = dq + 1/2 * dt * (v + v_new);
u_new(end/2+1:end) = v_new;


end
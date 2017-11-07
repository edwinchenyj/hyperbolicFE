function u_new = ImplicitMid( dt, u, obj, varargin)
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

it = 0;
dq_free = u(1:end/2);
v_free = u(end/2+1:end);
u_new = u;
v_free_new = u_new(end/2+1:end);


% initialize position vectors
% TODO: this should be move to inside the object
positionsM = obj.node';
positionsM = positionsM(:);
positions = positionsM;

positions(indLogical) = positionsM(indLogical) + dq_free + 1/4 * dt * (v_free + v_free_new);

% if isDG
    obj.SetCurrentState(positions - positionsM);
    K_mid = obj.StiffnessMatrix;
    Mass = obj.M;
    Eforce_mid = obj.ElasticForce;
% else
%     obj.SetCurrentState(positions - positionsM);
%     K_mid = obj.StiffnessMatrix;
%     Mass = obj.M;
%     Eforce_mid = obj.ElasticForce;
% end

Mass = Mass(indLogical,indLogical);
K_mid = K_mid(indLogical,indLogical);
B = -obj.a * Mass - obj.b * K_mid;

Eforce_mid = Eforce_mid(indLogical);

fExternal = Mass * obj.externalGravity;

f_mid = Eforce_mid + fExternal + B*1/2*(v_free + v_free_new);

N = obj.N;
nFixed = sum(constraint_indices);

residual0 = (dt * (Mass\f_mid))' * (dt * (Mass\f_mid));
Dv = -(speye(2*(N-nFixed)) + 1/4* dt*dt*(Mass\K_mid) - 1/2 * dt*(Mass\B))\(v_free_new - v_free - dt * (Mass\f_mid));
v_free_new = v_free_new + Dv;

residual = (v_free_new - v_free - dt * (Mass\f_mid))' * (v_free_new - v_free - dt * (Mass\f_mid));

it = it + 1;

u_new(1:end/2) = dq_free + 1/2 * dt * (v_free + v_free_new);
u_new(end/2+1:end) = v_free_new;

while (Dv'*Dv > 1e-12) && (residual > 1e-12)
    
    v_free_new = u_new(end/2+1:end);
    positions(indLogical) = positionsM(indLogical) + dq_free + 1/4*dt*(v_free + v_free_new);
    
%     if isDG
%         obj.SetCurrentDGState(positions - positionsM);
%     else
%         obj.SetCurrentState(positions - positionsM);
%     end
    
%     if isDG
        K_mid = obj.StiffnessMatrix;
        Mass = obj.M;
        Eforce_mid = obj.ElasticForce;
%     else
%         K_mid = obj.StiffnessMatrix;
%         Mass = obj.M;
%         Eforce_mid = obj.ElasticForce;
%     end
    
    Mass = Mass(indLogical,indLogical);
    K_mid = K_mid(indLogical,indLogical);
    
    B = -obj.a * Mass - obj.b * K_mid;
    
    Eforce_mid = Eforce_mid(indLogical);
    
    fExternal = Mass * obj.externalGravity;
    
    f_mid = Eforce_mid + fExternal + B*1/2*(v_free+v_free_new);
    
    Dv = -(speye(2*(N-nFixed)) + 1/4* dt*dt*(Mass\K_mid) - 1/2 * dt*(Mass\B))\(v_free_new - v_free - dt * (Mass\f_mid));
    v_free_new = v_free_new + Dv;
    
    residual = (v_free_new - v_free - dt * (Mass\f_mid))' * (v_free_new - v_free - dt * (Mass\f_mid));
    it = it + 1;
    
    u_new(1:end/2) = dq_free + 1/2 * dt * (v_free + v_free_new);
    u_new(end/2+1:end) = v_free_new;
    
    if (it > 3 && residual > residual0) || it == MaxIT
        disp('local substep required')
        if nargin > 3
            u_half = ImplicitMid(dt/2, u, obj, varargin);
        else
            u_half = ImplicitMid(dt/2, u, obj);
        end
        v_free = u_half(end/2 + 1:end);
        dq_free = u_half(1:end/2);
        
%         v(indLogical) = v_free;
        positions(indLogical) = positionsM(indLogical) + dq_free;
        
        u_half = [positions(indLogical)-positionsM(indLogical); v_free];
        
        if nargin > 3
            u_new = ImplicitMid(dt/2, u_half, obj, varargin);
        else
            u_new = ImplicitMid(dt/2, u_half, obj);
        end
        break;
    end
    
    if it == MaxIT
        error('Newton iteration not converging in IM')
    end
    
end

end
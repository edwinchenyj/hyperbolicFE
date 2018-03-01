function u_new = reduced_BackwardEuler( dt, u, obj, varargin)
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
obj.indLogical = indLogical;
Dim = obj.Dim;

it = 0;
dq = u(1:end/2);
v = u(end/2+1:end);
u_new = u;
v_new = u_new(end/2+1:end);

obj.x = obj.X + dq + dt * (v_new);

obj.SetCurrentState(obj.x - obj.X); % update deformation gradient
K_end = obj.StiffnessMatrix;
Mass = obj.M;
Mass = Mass(indLogical,indLogical);
K_end = K_end(indLogical,indLogical);



reduced_dimension = 100;
[V,D] = eigs(K_end,Mass,reduced_dimension,'sm');
[low_eig, permutation_indices] = sort(diag(D));
V = V(:,permutation_indices);


reduced_dq =  V'*Mass*(obj.x(indLogical)-obj.X(indLogical));
reduced_v =  V'*Mass*(v(indLogical));
reduced_v_new =  V'*Mass*(v_new(indLogical));

B = -obj.a * Mass - obj.b * K_end;

Eforce_end = obj.ElasticForce;
Eforce_end = Eforce_end(indLogical);

fExternal = Mass * obj.externalGravity(indLogical);

f_end = Eforce_end + fExternal + B*v_new(indLogical);

% f = Eforce + fExternal;

reduced_f =  V' * f_end;

reduced_residual0 = (dt * (reduced_f))' * (dt * (reduced_f));
reduced_Dv = -(speye(reduced_dimension) + dt*dt*(diag(low_eig)))\(reduced_v_new - reduced_v - dt * reduced_f);
% Dv = -(speye(reduced_dimension) + dt*dt*(diag(low_eig)) - dt*(Mass\B))\(v_new(indLogical) - v(indLogical) - dt * (Mass\f_end));
reduced_v_new = reduced_v_new + reduced_Dv;
reduced_dq = reduced_dq + dt * reduced_v_new;
% reduced_dq = redu

dq(indLogical) = V * reduced_dq; 
v(indLogical) = V * reduced_v;
u_new = [dq; v];

% v_new(indLogical) = v_new(indLogical) + Dv;

reduced_residual = (reduced_v_new - reduced_v - dt * reduced_f)' * (reduced_v_new - reduced_v - dt * reduced_f);

it = it + 1;

% u_new(1:end/2) = dq + dt * (v_new);
% u_new(end/2+1:end) = v_new;

% while (reduced_Dv'*reduced_Dv > 1e-12) && (reduced_residual > 1e-12)
%     
% %     v_new = u_new(end/2+1:end);
%     obj.x(indLogical) = obj.X(indLogical) + dq(indLogical) + dt * V * (reduced_v_new);
%     
%     obj.SetCurrentState(obj.x - obj.X);
%     K_end = obj.StiffnessMatrix;
% Mass = obj.M;
% Mass = Mass(indLogical,indLogical);
% K_end = K_end(indLogical,indLogical);
% 
% 
% 
% reduced_dimension = 10;
% [V,D] = eigs(K_end,Mass,reduced_dimension,'sm');
% [low_eig, permutation_indices] = sort(diag(D));
% V = V(:,permutation_indices);
% 
% 
% reduced_dq =  V'*Mass*(obj.x(indLogical)-obj.X(indLogical));
% reduced_v =  V'*Mass*(v(indLogical));
% reduced_v_new =  V'*Mass*(v_new(indLogical));
% 
% % B = -obj.a * Mass - obj.b * K_end;
% 
% Eforce_end = obj.ElasticForce;
% Eforce_end = Eforce_end(indLogical);
% 
% fExternal = Mass * obj.externalGravity(indLogical);
% 
% f_end = Eforce_end + fExternal + B*v_new(indLogical);
% 
% % f = Eforce + fExternal;
% 
% reduced_f =  V' * f_end;
% 
%     reduced_Dv = -(speye(reduced_dimension) + dt*dt*(diag(low_eig)))\(reduced_v_new - reduced_v - dt * reduced_f);
% % Dv = -(speye(reduced_dimension) + dt*dt*(diag(low_eig)) - dt*(Mass\B))\(v_new(indLogical) - v(indLogical) - dt * (Mass\f_end));
% reduced_v_new = reduced_v_new + reduced_Dv;
% reduced_dq = reduced_dq + dt * reduced_v_new;
% dq(indLogical) = V * reduced_dq; 
% v(indLogical) = V * reduced_v;
% u_new = [dq; v];
% 
% % v_new(indLogical) = v_new(indLogical) + Dv;
% 
% reduced_residual = (reduced_v_new - reduced_v - dt * reduced_f)' * (reduced_v_new - reduced_v - dt * reduced_f);

% it = it + 1;
%     
% %     u_new(1:end/2) = dq + dt * (v_new);
% %     u_new(end/2+1:end) = v_new;
%     
%     if (it > 3 && reduced_residual > reduced_residual0) || it == MaxIT
%         disp('local substep required')
%         if nargin > 3
%             u_half = reduced_BackwardEuler(dt/2, u, obj, varargin{1});
%         else
%             u_half = reduced_BackwardEuler(dt/2, u, obj);
%         end
%         v = u_half(end/2 + 1:end);
%         dq = u_half(1:end/2);
%         
%         obj.x = obj.X + dq;
%         
%         u_half = [obj.x-obj.X; v];
%         
%         if nargin > 3
%             u_new = reduced_BackwardEuler(dt/2, u_half, obj, varargin{1});
%         else
%             u_new = reduced_BackwardEuler(dt/2, u_half, obj);
%         end
%         break;
%     end
%     
%     if it == MaxIT
%         error('Newton iteration not converging in BE')
%     end
%     
% end

end
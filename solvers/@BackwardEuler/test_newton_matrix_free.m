function test_newton_matrix_free

%% Initialize the object
% a single tet
node = [ 100 0 0; 0 100 0; 0 0 100; 40 3 74];
elem = [1 2 3 4];
N = size(node,1);
Eobj = elasticTetObj(node, elem);

Y = 10; % Young's modululs
P = 0.4; % Poisson ratio
rho = 1; % density

Eobj.SetMaterial( Y, P, rho, 1, 1); % set the first tet to be neohookean
Eobj.finalize(); % finalize the material of the object

% set the initial state as completely undeformed state with zero energy
Dx = 5 * rand(3*N,1); % a rand displacement field
v = zeros(3*N,1);
Eobj.SetCurrentState(Dx,v);
u0 = [Dx;v]; % the initial state vector

%% The control sequence of the simulation
% TODO add more control for external forces
pullingForce = 10;
fExternal = zeros(size(Dx));
randDirection = rand(size(Dx));
fExternal = fExternal + pullingForce * randDirection;
g = @(t, u) [zeros(3*N,1); fExternal];
JV = @(t, v) JacobianVectorProduct(Eobj, v); % In general, Jacobian can be a function of time as well, so I put it here.
% But here the Jacobian is really only a function of u only


J = @(t, u) Jacobian(Eobj, u);
obj = BackwardEuler(u0', true, JV, g); % note: initial condition to the ode should be a row vector
obj.HJ = J;



h = 0.1;


for i = 1:10
    
obj.state = [obj.state; zeros(1,size(obj.state,2))];
t = obj.T;
obj.hT = [obj.hT h]; % Need this line here for obj.ImhA to access h later


% Newton's method
OldState = obj.state(obj.CurrentIndex,:)';
NewState = obj.state(obj.CurrentIndex,:)';

D1 = ones(size(obj.IC')); % initialize the Newton step vector to be something large
D2 = D1;
g = obj.Hg(t, OldState);


obj.KSMobj.SetHMV(@obj.ImhJ); % set the function to the matrix-vector product for the krylov subspace method

J = obj.HJ(t, OldState); % produce the Jacobian from the OldState
Eobj.SetCurrentState(OldState(1:3*N),OldState(3*N+1:end)); % set the state for the Eobj to OldState to produce the same Jacobian for the matrix-vector product

    % solve in matrix free form
    obj.KSMobj.SetWarmStart(NewState);
    [D1, flag] = obj.KSMobj.solve(-NewState+OldState+h*obj.HJV(t, NewState)+h*g);
    
    % solve with matrix inversion
    D2 = (obj.I - h*J)\(-NewState + OldState + h * J * NewState + h * g);
    
    if norm(D1 - D2) < 1e-4
        D = D1;
        NewState = NewState + D;
    else
        disp('D1 != D2')
    end
    
    obj.state(obj.CurrentIndex + 1,:) = NewState';
    
    obj.CurrentIndex = obj.CurrentIndex + 1;
    obj.T = obj.T + h;
end


end




function Jv = JacobianVectorProduct(Eobj, v)
N = Eobj.GetNNodes(); % number of nodes

    Jv = [v(3*N+1:end); -Eobj.ElasticForceDifferential(v(1:3*N))]

end

function J = Jacobian(Eobj, u)
    % Converting an ODE state to an elastic force of the Eobj
    N = Eobj.GetNNodes(); % number of nodes
    Eobj.SetCurrentState(u(1:3*N),u(3*N+1:end));
    J = [zeros(3*N), eye(3*N); -Eobj.StiffnessMatrix, zeros(3*N)]
end


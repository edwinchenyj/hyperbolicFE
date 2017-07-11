classdef ImplicitMidpoint < FirstOrderIVP
    %RK4 Solving a first-order IVP by implicit midpoint
    %   TODO: Detailed explanation goes here
    %   The implicit midpoint method solves the nonlinear equation
    %
    %  F(u_{n+1}) = $u_{n+1} - u_n - kf( \frac{u_n + u_{n+1}}{2}) $
    %
    % wher $ u_n$ is the current state of the system and $u_{n+1}$ is the
    % state of the system at the next time step (current time + k). In this
    % implementation the nonlinear equation is solved using Newton's method
    % with GMRES for the linear solve in the inner iteration.
    properties
        HJV; % function handle to the matrix vector product of Jacobian
        % used for matrix-free in the nonlinear solve if using Newton's method.
        % inputs: t, u (notice the
        % Jacobian is depending on the current configuration, if wish to
        % use modified newton's method by freezing the jacobain, need to
        % add additional parameter)
        HJou; % function handle to the Jacobian update routine that updates
        % the Jacobian implicitly by updating the state. That is, updating
        % J(t, u) in the matrix free version of J(t, u)*v. This is
        % necessary because we are not constructing the Jacobian explicitly
        % in the matrix free implementation
        % The input u will
        % update the state, and hence update the function pointed by HJV.
        HJ; % function handle to the Jacobian. inputs: t, u
        isMatrixFree = false;
        I; % the identity matrix that will be used in the Newton step
        MaxIT = 400
        NIT; % number of iterations required for each step
        KSMobj; % a KrylovSubspaceMethod object for matrix-free version
        isNewtonConverge = true;
    end
    
    methods
        
        function obj = ImplicitMidpoint(Hf, IC, varargin)
            % varargin list if matrix free:
            % (true, HJV, Hg, HJou)
            % varargin list if not matrix free:
            % (false, HJ, Hg)
            % also takes a boolean indicating if the nonlinear equation
            % should be solved by Newton's method in a matrix-free fashion.
            % If want to use the
            % matrix-free algorithm of the method, an additional function
            % handle to the matrix-vector product should be provided,
            % otherwise the function handle to the Jacobian of the function
            % should be provided. If none is provided, use other nonlinear
            % solvers like fixed-point, secant, or BFGS
            % use secant method by
            % default (NOT IMPLEMENTED YET!). TODO: implement secant method
            % TODO: BE doesn't require Hf, we only need HMV or HJ, and Hg!
            obj = obj@FirstOrderIVP(Hf, IC);
            

                
            obj.I = speye(size(IC,2));
            obj.state(1,:) = IC;
            obj.hT = [];
                
            if nargin > 2 % make sure varargin is given before using it
                assert(isa(varargin{2},'function_handle'));
                % varargin{1} is a boolean indicating if matrix-free or not
                if varargin{1}
                    % if isMatrixFree, the function handle is for the
                    % matrix-vector product
                    obj.HJV = varargin{2};
                    obj.isMatrixFree = true;
                    
                    % assert there are 5 inputs, with the last input as
                    % a function handle to the function to update the
                    % state (and hence the jacobian implicitly)
                    assert(nargin == 5);
                    assert(isa(varargin{3}, 'function_handle'));
                    obj.HJou = varargin{3};
                    
                    % the krylov subspace method that will be used for
                    % the matrix-free linear solve
                    obj.KSMobj = GMRES();
                else
                    % else the function handle is to the Jacobian
                    obj.HJ = varargin{2};
                    obj.isMatrixFree = false;
                end
                
            end
            
            
        end
        
        function sol = step(obj,k)
            %
            obj.state = [obj.state; zeros(1,size(obj.state,2))];
            t = obj.T;
            obj.hT = [obj.hT k]; % Need this line here for obj.ImhA to access k later
            
            
            % Newton's method
            OldState = obj.state(obj.CurrentIndex,:)'; % converting the row vector of the ODE state to a column vector for the linear solve 
            NewState = obj.state(obj.CurrentIndex,:)';
            
            D = zeros(size(obj.IC')); % initialize the Newton step
%             g = obj.Hg(t, OldState);
            it = 0;
            if obj.isMatrixFree
                % set the matrix for the krylove subspace matrix-free
                % linear solve to be (I - k/2*J(1/2(u_n + u_{n+1})), the
                % matrix for the implicit midpoint method
                obj.KSMobj.SetHMV(@obj.ImHalfkJ);
                
                % a "do-while" loop
                % update the Jacobian here only for inexact newton
                obj.HJou(t, (NewState+OldState)/2);
                % for the first iteration, (NewState+OldState)/2 = OldState                
                % set warm start for the krylov subspace method
                obj.KSMobj.SetWarmStart(D);
                
                D = obj.KSMobj.solve(-NewState+OldState + k * obj.Hf(obj.T,(NewState+OldState)/2));
                
                NewState = NewState + D;
                it = it + 1;
                while (D'*D > 1e-6) && (it < obj.MaxIT)
                    % update the Jacobian here for full newton step
                    obj.HJou(t, (NewState+OldState)/2);
                    
                    % set warm start for the krylov subspace method
                    obj.KSMobj.SetWarmStart(D);
                    
                    D = obj.KSMobj.solve(-NewState+OldState + k * obj.Hf(obj.T,(NewState+OldState)/2));
                    
                    NewState = NewState + D;
                    it = it + 1;
                    
                end
            else
                % inexact Newton: fixing the Jacobian and the nonlinear
                % part for each Newton iterations
                J = obj.HJ(t, (NewState+OldState)/2);
                % for the first iteration, (NewState+OldState)/2 = OldState
                D = (obj.I - 1/2*k*J)\(-NewState + OldState + k * obj.Hf(obj.T,(NewState+OldState)/2));
                NewState = NewState + D;
                it = it + 1;
                while (D'*D > 1e-6) && (it < obj.MaxIT)
                    % for full newton, update J in each iteration
                    J = obj.HJ(t, (NewState+OldState)/2);
                    D = (obj.I - 1/2*k*J)\(-NewState + OldState + k * obj.Hf(obj.T,(NewState+OldState)/2));
                    NewState = NewState + D;
                    it = it + 1;
                end
            end
            if it == obj.MaxIT
                obj.isNewtonConverge = false;
            end
            obj.NIT = [obj.NIT it];
            obj.state(obj.CurrentIndex + 1,:) = NewState';
            
            obj.CurrentIndex = obj.CurrentIndex + 1;
            obj.T = obj.T + k;
            sol = obj.state(obj.CurrentIndex,:);
        end
        
        function sol = solve(obj, T)
            % evolve the state along the time in T
            obj.T = T(1);
            for i = 1:length(T)-1
                k = T(i+1) - T(i);
                obj.step(k);
            end
            sol = obj.state;
        end
        
    end
    
    methods (Access = private)
        function MV = ImHalfkJ(obj,v)
            % matrix vector product for (I - 1/2*k*J)v. obj.HJV(v) gives the
            % matrix vector prodcut for J*v
            if size(v,1) == 1
                v = v'; % make sure v is a column vector
            end
            t = obj.T;
            MV = obj.I*v - 1/2*obj.hT(end) * obj.HJV(t, v);
        end
    end
    
    methods (Static)
        test_implicit_midpoint() % unit test on implicit midpoint
        test_implicit_midpoint_sine_gordon()
        test_implicit_midpoint_nonlinear_wave()
    end
    
end
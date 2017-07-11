classdef BackwardEuler < FirstOrderIVP
    %RK4 Solving a first-order IVP by backward Euler
    %   Detailed explanation goes here
    properties
        HJV; % function handle to the matrix vector product of Jacobian
        % used for matrix-free calculation. inputs: t, u (notice the
        % Jacobian is depending on the current configuration, if wish to
        % use modified newton's method by freezing the jacobain, need to
        % add additional parameter)
        HJou; % function handle to the Jacobian update routine that updates
        % the Jacobian implicitly by updating the state. The input u will
        % update the state, and hence update the function pointed by HJV.
        HJ; % function handle to the Jacobian. inputs: t, u
        isMatrixFree = false;
        I; % the identity matrix that will be used in the Newton step
        MaxIT = 40;
        KSMobj; % a KrylovSubspaceMethod object for matrix-free version
        isNewtonConverge = true;
    end
    
    methods
        
        function obj = BackwardEuler(Hf, IC, varargin)
            % varargin list if matrix free:
            % (true, HJV, Hg, HJou)
            % varargin list if not matrix free:
            % (false, HJ, Hg)
            % also takes a boolean indicating if the nonlinear equation
            % should be solved in a matrix-free fashion. If want to use the
            % matrix-free algorithm of the method, an additional function
            % handle to the matrix-vector product should be provided,
            % otherwise the function handle to the Jacobian of the function
            % should be provided. If none is provided, use secant method by
            % default (NOT IMPLEMENTED YET!). TODO: implement secant method
            % TODO: BE doesn't require Hf, we only need HMV or HJ, and Hg!
            obj = obj@FirstOrderIVP(Hf, IC);
            
            if nargin > 2 % make sure varargin is given
                
                obj.I = speye(size(IC,2));
                obj.state(1,:) = IC;
                obj.hT = [];
                
                
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
        
        function sol = step(obj,h)
            %
            obj.state = [obj.state; zeros(1,size(obj.state,2))];
            t = obj.T;
            obj.hT = [obj.hT h]; % Need this line here for obj.ImhA to access h later
            
            
            % Newton's method
            OldState = obj.state(obj.CurrentIndex,:)'; % converting the row vector of the ODE state to a column vector for the linear solve 
            NewState = obj.state(obj.CurrentIndex,:)';
            
            D = Inf(size(obj.IC')); % initialize the Newton step vector to be something large
%             g = obj.Hg(t, OldState);
            it = 0;
            if obj.isMatrixFree
                % set the matrix for the krylove subspace matrix-free
                % linear solve to be (I - h*J), the matrix for the backward
                % euler method
                obj.KSMobj.SetHMV(@obj.ImhJ);
                
                % update the Jacobian here for inexact newton
                % obj.HJou(t, OldState);
                while (D'*D > 1e-6) && (it < obj.MaxIT)
                    % update the Jacobian here for full newton step
                    obj.HJou(t, NewState);
                    
                    % set warm start for the krylov subspace method
                    if it ~= 0
                        obj.KSMobj.SetWarmStart(D);
                    else
                        obj.KSMobj.SetWarmStart(zeros(size(D)));
                    end
                    D = obj.KSMobj.solve(-NewState+OldState + h * obj.Hf(obj.T,NewState));
                    
                    NewState = NewState + D;
                    it = it + 1;
                    
                end
            else
                % inexact Newton: fixing the Jacobian and the nonlinear
                % part for each Newton iterations
                % J = obj.HJ(t, OldState);
                while (D'*D > 1e-6) && (it < obj.MaxIT)
                    % for full newton, update J in each iteration
                    J = obj.HJ(t, NewState);
                    D = (obj.I - h*J)\(-NewState + OldState + h * obj.Hf(obj.T,NewState));
                    NewState = NewState + D;
                    it = it + 1;
                end
            end
            if it == obj.MaxIT
                obj.isNewtonConverge = false;
            end
            
            obj.state(obj.CurrentIndex + 1,:) = NewState';
            
            obj.CurrentIndex = obj.CurrentIndex + 1;
            obj.T = obj.T + h;
            sol = obj.state(obj.CurrentIndex,:);
        end
        
        function sol = solve(obj, T)
            % evolve the state along the time in T
            obj.T = T(1);
            for i = 1:length(T)-1
                h = T(i+1) - T(i);
                obj.step(h);
            end
            sol = obj.state;
        end
        
    end
    
    methods (Access = private)
        function MV = ImhJ(obj,v)
            % matrix vector product for (I - h*J)v. obj.HMV(v) gives the
            % matrix vector prodcut for J*v
            if size(v,1) == 1
                v = v'; % make sure v is a column vector
            end
            t = obj.T;
            MV = obj.I*v - obj.hT(end) * obj.HJV(t, v);
        end
    end
    
    methods (Static)
%         test_newton_matrix_free() % unit test on newton's method, have access to all the properties, mehtods in the class
%         test_newton_matrix_free2()
        test_backward_euler() % unit test on backward euler
    end
    
end
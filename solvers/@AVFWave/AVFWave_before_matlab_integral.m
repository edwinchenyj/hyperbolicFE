classdef AVFWave < FirstOrderIVP
    % Solving a second-order IVP of wave equatino by AVF
    %   Detailed explanation goes here
    % TODO: make sure it's autonomous
    properties
        isMatrixFree = false;
        MaxIT = 400;
        KSMobj; % a KrylovSubspaceMethod object for matrix-free version
        H_intf; % the handle to the integrated RHS from the AVF method
        % main variable is delta_u, and the control parameters are u_n (u
        % at previous time step)
        NonlinearSolve = 0; % indicating which method to use for the nonlinear solve
        % 0: use matlab build-in fsolve
        % 1: use fixed-point iteration
        
        
        isUsingQuadrature = true;
        I; % identity of size of the system

        %% for newton's method
%         H_Dintf; % the handle to the Jacobian of the integrated RHS from the AVF method 
%         % i.e. Jacobian = H_Dintf(delta_u,u_n)
%         isNewtonConverge = true;
        
    end
    
    methods
        
        function obj = AVFWave(Hf, IC, NonlinearSolve, varargin)
            % If the integral of the RHS is not provided explicitly, use
            % quadrature rules to do numerical integration
            % varargin list if using quadrature and not matrix free:
            % {true}
            % varargin list if not using quadrature:
            % {false, H_intf}
            % varargin list if not using quadrature and using matrix free:
            % {false, H_intf, true, H_Dintf}
            % NonlinearSovle indicates which method to use for the nonlinear solve
            % 0: use matlab build-in fsolve
            % 1: use fixed-point iteration
            
            obj = obj@FirstOrderIVP(Hf, IC);
            
            obj.NonlinearSolve = NonlinearSolve;

            
            if nargin > 3
                if varargin{1}
                    obj.isUsingQuadrature = true;
                else
                    obj.isUsingQuadrature = false;
                    obj.H_intf = varargin{2};
                end
                
                if nargin > 5
                    if varargin{3}
                        obj.isMatrixFree = true;
%                         obj.H_Dintf = H_Dintf;
                    end
                end
            end
            obj.I = eye(length(IC));
        end
        
        function sol = step(obj,k)
            obj.state = [obj.state; zeros(1,size(obj.state,2))];
            obj.hT = [obj.hT k]; % the size of the time step
            t = obj.T;
            % currently using matlab fsolve to do the nonlinear solve
            % TODO: use Newton's method for AVF on wave
            % TODO: support matrix free-like equation
            u_n = obj.state(obj.CurrentIndex,:)';
            du0 = zeros(size(u_n));
            
            if obj.isUsingQuadrature
                quadrature_approximation = @(u_n, u_np1) 1/2 * (5/9 * obj.Hf(t,1/2*((1-sqrt(3/5))*u_n + (1+sqrt(3/5))*u_np1))...
                            + 8/9 * obj.Hf(t, 1/2 * (u_n + u_np1)) + 5/9 * obj.Hf(t,1/2*((1+sqrt(3/5))*u_n + (1-sqrt(3/5))*u_np1)));
                switch obj.NonlinearSolve
                    case 0
                        nonlinearEquation = @(delta_u) delta_u - k*quadrature_approximation(u_n,u_n+delta_u);
                        options = optimoptions('fsolve', 'TolFun', 1e-14, 'Display','off');
                        [du_sol, zeroNonlinearObjective, flag] = fsolve(nonlinearEquation, du0, options);
            
                    case 1
                        % fixed point iteration
                        du_sol = du0; % for checking with fixed point iteration
                        du_sol_old = du_sol;
                        fixedPointIteration = @(delta_u)  k*quadrature_approximation(u_n, u_n+delta_u);
                        du_sol = fixedPointIteration(du_sol);
                        it = 1;
                        while (norm(du_sol - du_sol_old) > 1e-9 && it < obj.MaxIT)
                            du_sol_old = du_sol;
                            du_sol = fixedPointIteration(du_sol);
                            it = it + 1;
                        end
                        if it >= obj.MaxIT
                            flag = 0;
                        else
                            flag = 1;
                        end
                    otherwise
                        % TODO: this is ugly... make logical switch
                        % prettier
                        quadrature_approximation = @(u_n, u_np1) 1/2 * (5/9 * obj.Hf(t,1/2*((1-sqrt(3/5))*u_n + (1+sqrt(3/5))*u_np1))...
                            + 8/9 * obj.Hf(t, 1/2 * (u_n + u_np1)) + 5/9 * obj.Hf(t,1/2*((1+sqrt(3/5))*u_n + (1-sqrt(3/5))*u_np1)));
                        nonlinearEquation = @(delta_u) delta_u - k*quadrature_approximation(u_n,u_n+delta_u);
                        options = optimoptions('fsolve', 'TolFun', 1e-14, 'Display','off');
                        [du_sol, zeroNonlinearObjective, flag] = fsolve(nonlinearEquation, du0, options);
                end
            else
                switch obj.NonlinearSolve
                    case 0 % use matlab build-in fsolve
                        nonlinearEquation = @(delta_u) delta_u - k*obj.H_intf(delta_u,u_n);
                        options = optimoptions('fsolve', 'TolFun', 1e-14, 'Display','off');
                        [du_sol, zeroNonlinearObjective, flag] = fsolve(nonlinearEquation, du0, options);
                        
                    case 1 % use fixed point iteration
                        
                        % fixed point iteration
                        du_sol = du0; % for checking with fixed point iteration
                        du_sol_old = du_sol;
                        fixedPointIteration = @(delta_u)  k*obj.H_intf(delta_u,u_n);
                        du_sol = fixedPointIteration(du_sol);
                        it = 1;
                        while (norm(du_sol - du_sol_old) > 1e-9 && it < obj.MaxIT)
                            du_sol_old = du_sol;
                            du_sol = fixedPointIteration(du_sol);
                            it = it + 1;
                        end
                        if it >= obj.MaxIT
                            flag = 0;
                        else
                            flag = 1;
                        end
                    otherwise
                        nonlinearEquation = @(delta_u) delta_u - k*obj.H_intf(delta_u,u_n);
                        options = optimoptions('fsolve', 'TolFun', 1e-14, 'Display','off');
                        [du_sol, zeroNonlinearObjective, flag] = fsolve(nonlinearEquation, du0, options);
                end
                
                % %                 trying newton
                %                 DnonlinearEquation = @(du) 1/k*obj.I - obj.H_Dintf(du,u);
                %                 newton = Newton(nonlinearEquation, delta_u, DnonlinearEquation, false);
                
                
            end
            
            if flag ~= 1
                warning('nonlinear solver not converging');
            end
            u_n = du_sol + u_n;
            
            obj.state(obj.CurrentIndex + 1,:) = u_n';
            
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
        
        function v_return = discreteGradient(obj, u_n, u_np1)
            t = obj.T;
            integrand = @(epsilon) obj.Hf(t, ( 1 - epsilon ) * u_n + epsilon * u_np1 );
            v_return = integral(integrand, 0, 1);
        end
    end

    methods (Static)
        test_AVFWave() % unit test on gautschi
    end

        
    
end
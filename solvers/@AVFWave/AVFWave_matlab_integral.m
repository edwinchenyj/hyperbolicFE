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
        
        
        Integration = 0; % indicate which integration method to use for the discrete gradient
        % 0: use 3 point Gauss Legendre quadrature
        % 1: use matlab build-in function integral
        % 2: use user provided analytical integrated function
        
        I; % identity of size of the system

        %% for newton's method
%         H_Dintf; % the handle to the Jacobian of the integrated RHS from the AVF method 
%         % i.e. Jacobian = H_Dintf(delta_u,u_n)
%         isNewtonConverge = true;
        
    end
    
    methods
        
        function obj = AVFWave(Hf, IC, NonlinearSolve, Integration, varargin)
            % If the integral of the RHS is not provided explicitly, use
            % Integration = 0 or 1 for numerical
            % integration. If the analytical integrated function is
            % provided, Integration = 2
            % varargin list if using numerical integral:
            % {} 
            % varargin list if numerical integral:
            % {H_intf}
            % NonlinearSovle indicates which method to use for the nonlinear solve
            % 0: use matlab build-in fsolve
            % 1: use fixed-point iteration
            % Integration indicates the integration method to use for the
            % discrete gradient        
            % 0: use 3 point Gauss Legendre quadrature
            % 1: use matlab build-in function integral
            % 2: use user provided analytical integrated function
            
            obj = obj@FirstOrderIVP(Hf, IC);
            if NonlinearSolve == 0 || NonlinearSolve == 1
                obj.NonlinearSolve = NonlinearSolve;
            else
                error('NonlinearSolve = 0 uses matlab build-in fsolve, 1 uses fixed-point iteration')
            end
            if Integration == 0 || Integration == 1 || Integration == 2
                obj.Integration = Integration;
            else
                error('Integration = 0 uses 3pt Gauss Legendre, 1 uses matlab build-in integral, 2 uses user provided integral')
            end
            
            if Integration == 2
                obj.H_intf = varargin{1};
            end
            obj.I = eye(length(IC));
        end
        
        function sol = step(obj,k)
            obj.state = [obj.state; zeros(1,size(obj.state,2))];
            obj.hT = [obj.hT k]; % the size of the time step
            t = obj.T;
            % TODO: use Newton's method for AVF on wave
            % TODO: support matrix free  equation
            u_n = obj.state(obj.CurrentIndex,:)';
            du0 = zeros(size(u_n));
            
            switch obj.NonlinearSolve
                case 0
                    nonlinearEquation = @(delta_u) delta_u - k * obj.discreteGradient(u_n, u_n + delta_u);
                    [du_sol, zeroNonlinearObjective, flag] = fsolve(nonlinearEquation, du0);
                case 1
                    % fixed point iteration
                    du_sol = du0; % for checking with fixed point iteration
                    du_sol_old = du_sol;
                    fixedPointIteration = @(delta_u)  k*obj.discreteGradient(u_n, u_n+delta_u);
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
            switch obj.Integration
                case 0
                    v_return = 1/2 * (5/9 * obj.Hf(t,1/2*((1-sqrt(3/5))*u_n + (1+sqrt(3/5))*u_np1))...
                        + 8/9 * obj.Hf(t, 1/2 * (u_n + u_np1)) + 5/9 * obj.Hf(t,1/2*((1+sqrt(3/5))*u_n + (1-sqrt(3/5))*u_np1)));
                case 1
                    integrand = @(epsilon) obj.Hf(t, ( 1 - epsilon ) * u_n + epsilon * u_np1 );
                    v_return = integral(integrand, 0, 1);
                case 2
                    v_return = obj.H_intf( u_n, u_np1 );
            end
        end
    end

    methods (Static)
        test_AVFWave() % unit test on gautschi
    end

        
    
end
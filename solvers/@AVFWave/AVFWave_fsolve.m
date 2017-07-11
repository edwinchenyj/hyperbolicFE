classdef AVFWave < FirstOrderIVP
    % Solving a second-order IVP of wave equatino by AVF
    %   Detailed explanation goes here
    % TODO: make sure it's autonomous
    properties
        isMatrixFree = false;
        MaxIT = 40;
        KSMobj; % a KrylovSubspaceMethod object for matrix-free version
        H_intf; % the handle to the integrated RHS from the AVF method
        % main variable is delta_u, and the control parameters are u_n (u
        % at previous time step)
        
        
        
        isUsingQuadrature = true;
        I; % identity of size of the system

        %% for newton's method
%         H_Dintf; % the handle to the Jacobian of the integrated RHS from the AVF method 
%         % i.e. Jacobian = H_Dintf(delta_u,u_n)
%         isNewtonConverge = true;
        
    end
    
    methods
        
        function obj = AVFWave(Hf, IC, varargin)
            % If the integral of the RHS is not provided explicitly, use
            % quadrature rules to do numerical integration
            % varargin list if using quadrature and not matrix free:
            % {true}
            % varargin list if not using quadrature:
            % {false, H_intf}
            % varargin list if not using quadrature and using matrix free:
            % {false, H_intf, true, H_Dintf}
            
            obj = obj@FirstOrderIVP(Hf, IC);
            if nargin > 2
                if varargin{1}
                    obj.isUsingQuadrature = true;
                else
                    obj.isUsingQuadrature = false;
                    obj.H_intf = varargin{2};
                end
                
                if nargin > 4
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
            obj.T = obj.T + k;
            t = obj.T;
            % currently using matlab fsolve to do the nonlinear solve
            % TODO: use Newton's method for AVF on wave
            % TODO: support matrix free-like equation
            u = obj.state(obj.CurrentIndex,:)';
            du0 = zeros(size(u));
            
            if obj.isUsingQuadrature
                quadrature_approximation = @(u_n, u_np1) 1/2 * (5/9 * obj.Hf(t,1/2*((1-sqrt(3/5))*u_n + (1+sqrt(3/5))*u_np1))...
                    + 8/9 * obj.Hf(t, 1/2 * (u_n + u_np1)) + 5/9 * obj.Hf(t,1/2*((1+sqrt(3/5))*u_n + (1-sqrt(3/5))*u_np1)));
                nonlinearEquation = @(delta_u) delta_u - k*quadrature_approximation(u,u+delta_u);
                options = optimoptions('fsolve', 'TolFun', 1e-14, 'Display','off');
                [du_sol, zeroNonlinearObjective, flag] = fsolve(nonlinearEquation, du0, options);
            else
                if obj.isMatrixFree
                    % TODO add matrix free version
                    
                else
                    nonlinearEquation = @(delta_u) delta_u - k*obj.H_intf(delta_u,u);
                    options = optimoptions('fsolve', 'TolFun', 1e-14, 'Display','off');
                    [du_sol, zeroNonlinearObjective, flag] = fsolve(nonlinearEquation, du0, options);
                    
                    % % Trying fixed point iteration
                    %                 du2 = du; % for checking with fixed point iteration
                    %             du2_old = du2;
                    %                 fixedPointIteration = @(delta_u)  k*obj.H_intf(delta_u,u);
                    %                 du2 = fixedPointIteration(du2);
                    %                 while (norm(du2 - du2_old) > 1e-6)
                    %                     du2_old = du2;
                    %                     du2 = fixedPointIteration(du2);
                    %                 end
                    %
                    %                 if norm(du2 - du) > 1e-4
                    %                     warning('fixed point not converging');
                    %                 end
                    
                    % %                 trying newton
                    %                 DnonlinearEquation = @(du) 1/k*obj.I - obj.H_Dintf(du,u);
                    %                 newton = Newton(nonlinearEquation, delta_u, DnonlinearEquation, false);
                end
                %             [flag, du] = newton.solve;
                
            end
            
            if flag ~= 1
                warning('nonlinear solver not converging');
            end
            u = du_sol + u;
            obj.state(obj.CurrentIndex + 1,:) = u';
            
            obj.CurrentIndex = obj.CurrentIndex + 1;
            
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

    methods (Static)
        test_AVFWave() % unit test on gautschi
    end

        
    
end
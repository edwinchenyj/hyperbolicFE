classdef ERAVF < FirstOrderIVP
    %ERAVF Solving a first-order IVP by ERAVF
    %(Rosenbrock-Euler)
    %   Detailed explanation goes here
    properties
        HJ; % function handle to the Jacobian. inputs: t, u
        Jn; % the Jacobian at each time step
        KSM; % krylov subspace method for matrix exponential approximation
        % TODO: to support true matrix free, need the matrix-vector product
        % with the jacobian
        SubspaceDimension = 15;
        isUsingKSM;
        
        NonlinearSolve = 0; % indicating which method to use for the nonlinear solve
        % 0: use matlab build-in fsolve
        % 1: use fixed-point iteration
        MaxIT = 40; % for fixed point approximation (if fixed point approximation is used)
        NIT; % number of fiexd point iteration for each step
        
        Integration = 0; % indicate which integration method to use for the discrete gradient
        % 0: use 3 point Gauss Legendre quadrature
        % 1: use matlab build-in function integral
    end

    methods
        function obj = ERAVF(Hf, IC, HJ, NonlinearSolve, Integration, usingKSM)
            if nargin == 0
                error('Please provide the function handle to the derivative and the initial condition')
            end
            obj = obj@FirstOrderIVP(Hf, IC);
            obj.HJ = HJ;
            obj.NonlinearSolve = NonlinearSolve;
            obj.Integration = Integration;
            if usingKSM
                obj.isUsingKSM = true;
            else
                obj.isUsingKSM = false;
            end
        end
        
        function sol = step(obj,k)
            obj.state = [obj.state; zeros(1,size(obj.state,2))]; % increase the size of obj.state
            
            t = obj.T;
            obj.hT = [obj.hT; k];
            u_n = obj.state(obj.CurrentIndex,:)'; % get a column vector of the current state
            du0 = zeros(size(u_n));
            % calculate the linear part if using KSM approx.
            % calculate exp(k*J)*u_n
            
            % temporarily store the jacobian at this time step
            obj.Jn = obj.HJ(t,u_n);
            
            if obj.isUsingKSM
                MV = @(v) obj.Jn * v;
                v0 = u_n;
                obj.KSM = Arnoldi(MV, v0, obj.SubspaceDimension);
                [Vmp1, HmBar, D] = obj.KSM.ConstructBasis;
                if D == (obj.SubspaceDimension)+1
                    Hm = HmBar(1:end-1,:);
                    Vm = Vmp1(:,1:end-1);
                else
                    Hm = HmBar(1:D,1:D);
                    Vm = Vmp1(:,1:D);
                end
                projected_exphJ = expm(k*Hm);
                E1 = zeros(size(projected_exphJ,2),1); E1(1) = 1;
                linear_part = norm(v0) * Vm * projected_exphJ * E1;
            else
                % skip for now
%                 linear_part = expm(k*obj.Jn) * u_n;
            end
            
            switch obj.NonlinearSolve
                case 0
                    nonlinearEquation = @(delta_u) delta_u + u_n - linear_part - k * obj.nonlinearDiscreteGradient(u_n, u_n + delta_u);
                    options = optimoptions('fsolve', 'TolFun', 1e-14, 'Display','off');
                    [du_sol, zeroNonlinearObjective, flag] = fsolve(nonlinearEquation, du0, options);
                case 1
                    % fixed point iteration
                    du_sol = du0; % for checking with fixed point iteration
                    du_sol_old = du_sol;
                    fixedPointIteration = @(delta_u)  -u_n + linear_part + k*obj.nonlinearDiscreteGradient(u_n, u_n+delta_u);
                    du_sol = fixedPointIteration(du_sol);
                    it = 1;
                    while (norm(du_sol - du_sol_old) > 1e-6 && it < obj.MaxIT)
                        du_sol_old = du_sol;
                        du_sol = fixedPointIteration(du_sol);
                        it = it + 1;
                    end
                    if it >= obj.MaxIT
                        flag = 0;
                    else
                        flag = 1;
                        obj.NIT = [obj.NIT it];
                    end
                case 2
                    m_max = [55];
                    p_max = [8];
                    u_np1 = u_n;
%                     u_np1(1:end/2) = u_np1(1:end/2) + k*u_n(end/2 + 1:end);
                    gn_tilde = (1/2 * (5/9 * obj.g( t , 1/2*((1-sqrt(3/5))*u_n + (1+sqrt(3/5))*u_np1))...
                        + 8/9 * obj.g( t , 1/2 * (u_n + u_np1)) + ...
                        5/9 * obj.g( t , 1/2*((1+sqrt(3/5))*u_n + (1-sqrt(3/5))*u_np1))));
                    A_tilde = sparse(size(obj.Jn,1) + 1, size(obj.Jn,1) + 1);
                    A_tilde(1:end-1,:) = [obj.Jn, gn_tilde];
                    u_tilde = [u_n; 1];
                    
                    M = select_taylor_degree(A_tilde,u_tilde,m_max,p_max);
        
        
                    [X, s, m] = expmv(k, A_tilde, u_tilde, M);
                    u_np1 = X(1:end-1);
                    du_sol = u_np1 - u_n;
                    du_sol_old = zeros(size(du_sol));
                    it = 1;
                    while (norm(du_sol - du_sol_old) > 1e-6 && it < obj.MaxIT)
                        du_sol_old = du_sol;
                        gn_tilde = (1/2 * (5/9 * obj.g( t , 1/2*((1-sqrt(3/5))*u_n + (1+sqrt(3/5))*u_np1))...
                        + 8/9 * obj.g( t , 1/2 * (u_n + u_np1)) + ...
                        5/9 * obj.g( t , 1/2*((1+sqrt(3/5))*u_n + (1-sqrt(3/5))*u_np1))));
                        A_tilde(1:end-1,end) = gn_tilde;
%                         M = select_taylor_degree(A_tilde,u_tilde,m_max,p_max);
                        [X, s, m] = expmv(k, A_tilde, u_tilde, M);
                        u_np1 = X(1:end-1);
                        du_sol = u_np1 - u_n;
                        it = it + 1;
                    end
                    if it >= obj.MaxIT
                        flag = 0;
                    else
                        flag = 1;
                        obj.NIT = [obj.NIT it];
                    end
                case 3
                    u_np1 = u_n;
                    u_np1(1:end/2) = u_np1(1:end/2) + k*u_n(end/2 + 1:end);
                    gn_tilde = (1/2 * (5/9 * obj.g( t , 1/2*((1-sqrt(3/5))*u_n + (1+sqrt(3/5))*u_np1))...
                        + 8/9 * obj.g( t , 1/2 * (u_n + u_np1)) + ...
                        5/9 * obj.g( t , 1/2*((1+sqrt(3/5))*u_n + (1-sqrt(3/5))*u_np1))));
                    A_tilde = sparse(size(obj.Jn,1) + 1, size(obj.Jn,1) + 1);
                    A_tilde(1:end-1,:) = [obj.Jn, gn_tilde];
                    u_tilde = [u_n; 1];
                    
        
                    [X] = expv(k, A_tilde, u_tilde);
                    u_np1 = X(1:end-1);
                    du_sol = u_np1 - u_n;
                    du_sol_old = zeros(size(du_sol));
                    it = 1;
                    while (norm(du_sol - du_sol_old) > 1e-6 && it < obj.MaxIT)
                        du_sol_old = du_sol;
                        gn_tilde = (1/2 * (5/9 * obj.g( t , 1/2*((1-sqrt(3/5))*u_n + (1+sqrt(3/5))*u_np1))...
                        + 8/9 * obj.g( t , 1/2 * (u_n + u_np1)) + ...
                        5/9 * obj.g( t , 1/2*((1+sqrt(3/5))*u_n + (1-sqrt(3/5))*u_np1))));
                        A_tilde(1:end-1,end) = gn_tilde;
                        [X, s, m] = expv(k, A_tilde, u_tilde);
                        u_np1 = X(1:end-1);
                        du_sol = u_np1 - u_n;
                        it = it + 1;
                    end
                    if it >= obj.MaxIT
                        flag = 0;
                    else
                        flag = 1;
                        obj.NIT = [obj.NIT it];
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
    end
    
    methods (Static)
        test_ERAVF_sine_gordon()
        test_ERAVF_sine_gordon_2()
        test_ERAVF_sine_gordon_3()
        test_ERAVF_nonlinear_wave()
    end
    
    methods (Access = protected)
        
        function A = phi1(obj,z)
%             [V,D] = eig(full(z));
%             for i = 1:length(z)
%                     D(i,i) = (expm1(D(i,i)))/D(i,i);
%             end
%             A = V*D/(V);
%             A = (expm(z) - eye(size(z)))/z;
                        [V,D] = eig(full(z));
            for i = 1:length(z)
                    if (D(i,i) ~=0)
                        D(i,i) = (expm1(D(i,i)))/D(i,i);
                    else
                        D(i,i) = 1;
                    end
            end
            A = V*D/(V);
        end

        function v_return = g(obj, t, v)
            v_return = obj.Hf(t, v) - obj.Jn * v;
        end
        
        function v_return = nonlinearDiscreteGradient(obj, u_n, u_np1)
            % discrete gradient of the nonlinear part
            t = obj.T;
            k = obj.hT(end);
            switch obj.Integration
                case 0
                    v_return = real(obj.phi1( k * obj.Jn )) * ...
                        (1/2 * (5/9 * obj.g( t , 1/2*((1-sqrt(3/5))*u_n + (1+sqrt(3/5))*u_np1))...
                        + 8/9 * obj.g( t , 1/2 * (u_n + u_np1)) + ...
                        5/9 * obj.g( t , 1/2*((1+sqrt(3/5))*u_n + (1-sqrt(3/5))*u_np1))));
                case 1
                    integrand = @(epsilon) obj.g( t, ( 1 - epsilon ) .* u_n + epsilon .* u_np1 );
                    v_return = real(obj.phi1( k * obj.Jn )) * ...
                        integral(integrand, 0, 1, 'ArrayValued', true);
            end
        end
    end
end


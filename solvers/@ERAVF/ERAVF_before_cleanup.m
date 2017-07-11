classdef ERAVF < FirstOrderIVP
    %ERAVF Solving a first-order IVP by ERAVF
    %(Rosenbrock-Euler)
    %   Detailed explanation goes here
    properties
        HJ; % function handle to the Jacobian. inputs: t, u
        KSM; % krylov subspace method for matrix exponential approximation
        % TODO: to support true matrix free, need the matrix-vector product
        % with the jacobian
        SubspaceDimension = 15;
        isUsingKSM;
    end

    methods
        function obj = ERAVF(Hf, IC, HJ, usingKSM)
            if nargin == 0
                error('Please provide the function handle to the derivative and the initial condition')
            end
            obj = obj@FirstOrderIVP(Hf, IC);
            if nargin >= 3
                obj.HJ = HJ;
                obj.isUsingKSM = usingKSM;
            end
        end
        
        function sol = step(obj,k)
            assert(isa(obj.Hf,'function_handle'));
            obj.state = [obj.state; zeros(1,size(obj.state,2))]; % increase the size of obj.state
            obj.T = obj.T + k;
            t = obj.T;
            obj.hT = [obj.hT; k];
            u_n = obj.state(obj.CurrentIndex,:)'; % get a column vector of the current state
            MV = @(v) obj.HJ(t,u_n) * v;
            du0 = zeros(size(u_n));
            % calculate the linear part if using KSM approx.
            % calculate exp(k*J)*u_n
            if obj.isUsingKSM
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
                projected_exphJ = obj.matrixExp(k*Hm);
                E1 = zeros(size(projected_exphJ,2),1); E1(1) = 1;
                linear_part = norm(v0) * Vm * projected_exphJ * E1;
                nonlinearEquation = @(delta_u) delta_u + u_n ...
                    - linear_part - k * obj.NonlinearPartQuadrature(u_n, u_n + delta_u);
            else
                linear_part = expm(k*obj.HJ(t,u_n)) * u_n;
                nonlinearEquation = @(delta_u) delta_u + u_n ...
                    - linear_part - k * obj.NonlinearPartQuadrature(u_n, u_n + delta_u);
            end
            options = optimoptions('fsolve', 'TolFun', 1e-9, 'Display','final-detailed');
            [du_sol, zeroNonlinearObjective, flag] = fsolve(nonlinearEquation, du0, options);
            
            if flag ~= 1
                warning('nonlinear solver not converging');
            end
            new_u = du_sol + u_n;
            
            obj.state(obj.CurrentIndex + 1,:) = new_u';
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
        test_ERAVF()
    end
    
    methods (Access = protected)
        function A = matrixExp(obj,z)
            [V,D] = eig(full(z));

            A = V*diag(exp(diag(D)))/V;
        end
        
        function A = phi1(obj,z)
            [V,D] = eig(full(z));
            for i = 1:length(z)
                if D(i,i) ~= 0
                    D(i,i) = (exp(D(i,i)) - 1)/D(i,i);
                else
                    D(i,i) = 1;
                end
            end
            A = V*D/V;
        end

        function v_return = g(obj, t, v)
            v_return = obj.Hf(t, v) - obj.HJ(t, v) * v;
        end
        
        function v_return = NonlinearPartQuadrature_KSM(obj, u_n, u_np1)
            
            % TODO: support true matrix free calculation when using
            % calculate the nonlinear part using AVF
            t = obj.T;
            MV = @(v) obj.HJ(t,u_n) * v;

            v0 = 1/2 * (5/9 * obj.g(t,1/2*((1-sqrt(3/5))*u_n + (1+sqrt(3/5))*u_np1))...
                + 8/9 * obj.g(t, 1/2 * (u_n + u_np1)) + 5/9 * obj.g(t,1/2*((1+sqrt(3/5))*u_n + (1-sqrt(3/5))*u_np1)));
            
            obj.KSM = Arnoldi(MV, v0, obj.SubspaceDimension);
            [Vmp1, HmBar, D] = obj.KSM.ConstructBasis;
            if D == (obj.SubspaceDimension)+1
                Hm = HmBar(1:end-1,:);
                Vm = Vmp1(:,1:end-1);
            else
                Hm = HmBar(1:D,1:D);
                Vm = Vmp1(:,1:D);
            end
            k = obj.hT(end);
            %%%%% should do eig first eig(Hm)
            projected_phi1 = obj.phi1(k*Hm);
            E1 = zeros(size(projected_phi1,2),1); E1(1) = 1;
            v_return = norm(v0) * Vm * projected_phi1 * E1;
            
            
        end
        
        function v_return = NonlinearPartQuadrature(obj, u_n, u_np1)
            k = obj.hT(end);
            t = obj.T;
            v_return = obj.phi1(k*obj.HJ(t,u_n)) * 1/2 * (5/9 * obj.g(t, 1/2 * ( ( 1-sqrt(3/5) ) * u_n + (   1+sqrt(3/5) )*u_np1 ) )...
                + 8/9 * obj.g(t, 1/2 * ( u_n + u_np1) )...
                + 5/9 * obj.g(t, 1/2 * ( ( 1+sqrt(3/5) ) * u_n + ( 1-sqrt(3/5) ) * u_np1 ) ) );
            
        end
    end
end


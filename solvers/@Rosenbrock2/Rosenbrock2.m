classdef Rosenbrock2 < FirstOrderIVP
    %Rosenbrock2 Solving a first-order IVP by 2nd order Rosenbrock
    %(Rosenbrock-Euler)
    %   Detailed explanation goes here
    properties
        HJ; % function handle to the Jacobian. inputs: t, u
        Jn; % Jacobian for each time step
        KSM; % krylov subspace method for matrix exponential approximation
        % TODO: to support true matrix free, need the matrix-vector product
        % with the jacobian
        SubspaceDimension = 40;
        isUsingKSM;
    end
    
    methods
        function obj = Rosenbrock2(Hf, IC, HJ, usingKSM)
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
            u_n = obj.state(obj.CurrentIndex,:)'; % get a column vector of the current state
            
            obj.Jn = obj.HJ(t, u_n);
            g = obj.g(t,u_n);
            eta = 2 ^ (-ceil(log2(norm(g,1))));
            norm(eta*g,1)
            A_tilde = sparse(size(obj.Jn,1) + 1, size(obj.Jn,1) + 1);
            A_tilde(1:end-1,:) = [obj.Jn, eta*g];
            u_tilde = [u_n; 1/eta];
            X = expmv(k, A_tilde,u_tilde);
            if isnumeric(X)
                u_n1 = X(1:end-1);
                
                obj.state(obj.CurrentIndex + 1,:) = u_n1;
                %             if obj.isUsingKSM
                %                 % TODO: support true matrix free calculation when using
                %                 % Krylov subspace method to do matrix approximation
                %                 % calculate exp(k*J)*u_n
                %                 MV = @(v) obj.Jn * v;
                %                 v0 = u_n;
                %                 obj.KSM = Arnoldi(MV, v0, obj.SubspaceDimension);
                %                 [Vmp1, HmBar, D] = obj.KSM.ConstructBasis;
                %                 if D == (obj.SubspaceDimension)+1
                %                     Hm = HmBar(1:end-1,:);
                %                     Vm = Vmp1(:,1:end-1);
                %                 else
                %                     Hm = HmBar(1:D,1:D);
                %                     Vm = Vmp1(:,1:D);
                %                 end
                %                 projected_exphJ = obj.matrixExp(k*Hm);
                %                 E1 = zeros(size(projected_exphJ,2),1); E1(1) = 1;
                %                 linear_part = norm(v0) * Vm * projected_exphJ * E1;
                %                 % calculate phi1(k*J) * g(t,u) , the nonlinear part
                %                 MV = @(v) obj.Jn * v;
                %                 v0 = obj.g(t,u_n);
                %                 if norm(v0) < 1e-12
                %                     obj.state(obj.CurrentIndex + 1,:) =  (linear_part)';
                %                 else
                %                     obj.KSM = Arnoldi(MV, v0, obj.SubspaceDimension);
                %                     [Vmp1, HmBar, D] = obj.KSM.ConstructBasis;
                %                     if D == (obj.SubspaceDimension)+1
                %                         Hm = HmBar(1:end-1,:);
                %                         Vm = Vmp1(:,1:end-1);
                %                     else
                %                         Hm = HmBar(1:D,1:D);
                %                         Vm = Vmp1(:,1:D);
                %                     end
                %                     projected_phi1 = real(obj.phi1(k*Hm));
                %                     E1 = zeros(size(projected_phi1,2),1); E1(1) = 1;
                %                     obj.state(obj.CurrentIndex + 1,:) =  (linear_part + k * norm(v0) * Vm * projected_phi1 * E1)';
                %                 end
                %             else
                %                 obj.state(obj.CurrentIndex + 1,:) =  (u_n + k * real(obj.phi1(k*obj.Jn))*obj.Hf(t,u_n))';
                %             end
                
                
                obj.CurrentIndex = obj.CurrentIndex + 1;
                obj.hT = [obj.hT; k];
                sol = obj.state(obj.CurrentIndex,:);
            else
                sol = X; % output the flag from expmv
            end
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
        test_rosenbrock2()
        test_rosenbrock2_sine_gordon()
    end
    
    methods
        
        function A = matrixExp(obj,z)
            [V,D] = eig(full(z));
            
            A = V*diag(exp(diag(D)))/(V);
        end
        function A = phi1(obj,z)
            [V,D] = eig(full(z));
            for i = 1:length(z)
                D(i,i) = (expm1(D(i,i)))/D(i,i);
            end
            A = V*D/(V);
        end
        
        function v_return = g(obj, t, v)
            v_return = obj.Hf(t, v) - obj.Jn * v;
        end
    end
end


classdef Rosenbrock2 < FirstOrderIVP
    %Rosenbrock2 Solving a first-order IVP by 2nd order Rosenbrock
    %(Rosenbrock-Euler)
    %   Detailed explanation goes here
    properties
        HJ; % function handle to the Jacobian. inputs: t, u
        KSM; % krylov subspace method for matrix exponential approximation
        % TODO: to support true matrix free, need the matrix-vector product
        % with the jacobian
        SubspaceDimension = 15;
    end

    methods
        function obj = Rosenbrock2(Hf, IC, HJ)
            if nargin == 0
                error('Please provide the function handle to the derivative and the initial condition')
            end
            obj = obj@FirstOrderIVP(Hf, IC);
            if nargin == 3
                obj.HJ = HJ;
            end
        end
        
        function sol = step(obj,h)
            assert(isa(obj.Hf,'function_handle'));
            obj.state = [obj.state; zeros(1,size(obj.state,2))]; % increase the size of obj.state
            t = obj.T;
            st = obj.state(obj.CurrentIndex,:)'; % get a column vector of the current state
            % TODO: support true matrix free multiplication
            % calculate exp(k*J)*st
            MV = @(v) obj.HJ(t,st) * v;
            v0 = st;
            obj.KSM = Arnoldi(MV, v0, obj.SubspaceDimension);
            [Vmp1, HmBar, D] = obj.KSM.ConstructBasis;
            if D == (obj.SubspaceDimension)+1
                Hm = HmBar(1:end-1,:);
                Vm = Vmp1(:,1:end-1);
            else
                Hm = HmBar(1:D,1:D);
                Vm = Vmp1(:,1:D);
            end
            projected_exphJ = expm(h*Hm);
            E1 = zeros(size(projected_exphJ,2),1); E1(1) = 1;
            linear_part = norm(v0) * Vm * projected_exphJ * E1;
            % calculate phi1(k*J) * g(t,u) , the nonlinear part
            MV = @(v) obj.HJ(t,st) * v;
            v0 = obj.g(t,st);
            obj.KSM = Arnoldi(MV, v0, obj.SubspaceDimension);
            [Vmp1, HmBar, D] = obj.KSM.ConstructBasis;
            if D == (obj.SubspaceDimension)+1
                Hm = HmBar(1:end-1,:);
                Vm = Vmp1(:,1:end-1);
            else
                Hm = HmBar(1:D,1:D);
                Vm = Vmp1(:,1:D);
            end
            projected_phi1 = obj.phi1(h*Hm);
            E1 = zeros(size(projected_phi1,2),1); E1(1) = 1;
            obj.state(obj.CurrentIndex + 1,:) =  (linear_part + h * norm(v0) * Vm * projected_phi1 * E1)';
%             obj.state(obj.CurrentIndex + 1,:) =  (st + h * obj.phi1(h*obj.HJ(t,st))*obj.Hf(t,st))';
%             obj.state(obj.CurrentIndex + 1,:) =  (st + h * obj.phi1(h*obj.HJ(t,st))* obj.Hf(t,st))';
            
            obj.CurrentIndex = obj.CurrentIndex + 1;
            obj.T = obj.T + h;
            obj.hT = [obj.hT; h];
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
    
    methods (Static)
        test_rosenbrock2()
    end
    
    methods (Access = protected)
        function A = phi1(obj,z)
            [V,D] = eig(full(z));
            for i = 1:length(z)
                if (D(i,i) ~=0)
                    D(i,i) = (exp(D(i,i)) - 1)/D(i,i);
                else
                    D(i,i) = 1;
                end
            end
            A = V*D*inv(V);
        end
        
        function A = phi2(obj,z)
            A = z\(phi1(z) - eye(size(z)));
        end
        function v_return = g(obj, t, v)
            v_return = obj.Hf(t, v) - obj.HJ(t, v) * v;
        end
    end
end


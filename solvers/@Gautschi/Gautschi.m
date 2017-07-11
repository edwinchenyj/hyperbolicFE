classdef Gautschi < SecondOrderIVP
    %RK4 Solving a first-order IVP by backward Euler
    %   Detailed explanation goes here
    % TODO: make sure it's autonomous
    properties
        isMatrixFree = false;
        MaxIT = 40;
        KSMobj; % a KrylovSubspaceMethod object for matrix-free version
        isNewtonConverge = true;
        Omega;
        Omega2; % Omega^2 from Hochbrock and Oosterman which 
        Argumentfiltering = false; % the option for nonlinear system
        % is equal to the A in Michels
        
        isAxFactorized = false;
        isFilteredAxFactorized = false;
        
        phiFilter = true; % true: phi = sinc, false: phi = 1
        psiFilter = true; % true: psi = sinc * phi, false: psi = phi 
        % two flags used to save computation since some Arnoldi iterations
        % and eigenvalue decomposition are operating on the same matrix
        % vector pair
        
        O2_V; % [O2_V, O2_D] = eig(Omega2)
        O2_D;
        isSinCal = false; % flags to indicate whether the matrix products are calculated or not
        isOsinCal = false; % this is omega * sin(t*omega)
        isCosCal = false;
        isSincCal = false;
        isSinc2Cal = false;
        SinCal; % variables to store calculated matrices
        OsinCal;
        CosCal
        SincCal;
        Sinc2Cal;
    end
    
    methods
        
        function obj = Gautschi(HD, HDJ, IC, M, varargin)
            % varargin list if matrix free:
            % {HB, true}
            % varargin list if not matrix free:
            % {HB, false}
            
            obj = obj@SecondOrderIVP(HD, HDJ, IC, M, varargin);
            if nargin == 6
                if varargin{2}
                    obj.isMatrixFree = true;
                end
            end
        end
        
        function sol = step(obj,h)
            obj.isOsinCal = false;
            obj.isSinCal = false; % set the flags down for each new state
            obj.isCosCal = false;
            obj.isSincCal = false;
            obj.isSinc2Cal = false;
            
            obj.state = [obj.state; zeros(1,size(obj.state,2))];
            t = obj.T;
            obj.hT = [obj.hT h]; % the size of the time step
            
            q = obj.state(obj.CurrentIndex,1:(size(obj.state,2)/2))';
            p = obj.state(obj.CurrentIndex,(size(obj.state,2)/2)+1:end)';
            % TODO: support matrix free
            obj.Omega2 =  obj.HOmega2(t,q);% Omega^2

            
            if length(obj.hT) ~= 1
                q_old = obj.state(obj.CurrentIndex-1,1:(size(obj.state,2)/2))';
                % if it's not the first step, use the 2-step method
                q_phi = MatrixFuncVectorProduct(0, obj.Omega2, q);
                    if obj.Argumnentfiltering
                        obj.Omega2 = obj.HOmega2(t,q_phi);
                        q_phi = MatrixFuncVectorProduct(0, obj.Omega2, q); % 0 corresponds to phi
                    end
                v_phi = 1/t * (q - q_old);
                nonlinear_part_phi = obj.Mg(q_phi);
                q_c = MatrixFuncVectorProduct(1, obj.Omega2, q); % 1 corresponds to cos
                q_psi = MatrixFuncVectorProduct(2, obj.Omega2, nonlinear_part_phi); % 2 corresponds to psi
                new_q = 2 * q_c - q_old - obj.hT(end)^2 * obj.M\q_psi;        
                new_p = 1/t * (new_q - q);
            else 
                % if it's the first step, use the 1-step method
                if obj.isMatrixFree
                    % TODO add matrix free version
                    
                else
                    % R = [obj.cos_hTOmega, obj.sinc_hTOmega...
                    %      -obj.osin_hTOmega, obj.cos_hTOmega];
                    [obj.O2_V, obj.O2_D] = eig(obj.Omega2);
                    new_q = [obj.cos_hTOmega, h*obj.sinc_hTOmega] * [q;p] + h/2 * h * obj.psi_hTOmega * obj.g(obj.phi_hTOmega * q);
                    new_p = [-obj.osin_hTOmega, obj.cos_hTOmega] * [q;p] + h/2 * obj.psi0_hTOmega * obj.g(obj.phi_hTOmega * q) + obj.psi1_hTOmega * obj.g(obj.phi_hTOmega * new_q);
                end
            

            end
            obj.state(obj.CurrentIndex + 1,:) = [new_q; new_p]';
            
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
    
    
    
    methods (Static)
        test_gautschi() % unit test on gautschi

    end

    methods (Access = protected)
        function f = phi(obj,x)
            if obj.phiFilter
                f = sinc(x/pi);
            else
                f = x;
            end
        end
        function f = psi(obj,x)
            if obj.psiFilter
                f = sinc(x/pi) * phi(x);
            else
                f = phi(x);
            end
        end
        
        function v_return = MatrixFuncVectorProduct(obj,FuncType, A, v)
            MV = @(v) A*v;
            m = obj.SubspaceDimension;
            arnoldi = Arnoldi(MV,v,m);
            [Vmp1, HmBar, D] = arnoldi.ConstructBasis;
            if D == m+1
                Hm = HmBar(1:end-1,:);
                Vm = Vmp1(:,1:end-1);
            else
                Hm = HmBar(1:D,1:D);
                Vm = Vmp1(:,1:D);
            end
            [Evec, DiagEval] = eig(Hm);
            switch FuncType

                case 0 % case 0 corresponds to phi(k * sqrt(A))*v
                    for i = 1:m
                        DiagEval(i,i) = obj.phi(obj.hT(end) * sqrt(DiagEval(i,i)));
                    end
                case 1 % case 1 corresponds to cos(k * sqrt(A))*v
                    for i = 1:m
                        DiagEval(i,i) = cos(obj.hT(end) * sqrt(DiagEval(i,i)));
                    end
                case 2 % case 2 corresponds to psi(k * sqrt(A))*v
                    for i = 1:m
                        DiagEval(i,i) = obj.psi(obj.hT(end) * sqrt(DiagEval(i,i)));
                    end
            end
            e1 = zeros(m,1); e1(1) = 1;
            z = Evec\DiagEval*Evec * e1;
            v_return = norm(v) * Vm * z;
        end
        
        function A = cos_hTOmega(obj)
            if ~obj.isCosCal
                z_v = obj.O2_V;
                z_d = obj.O2_D;
                h = obj.hT(end);
                for i = 1:length(z_v)
                    if(z_d(i,i) ~= 0)
                        z_d(i,i) = cos(h*sqrt(z_d(i,i)));
                    else
                        z_d(i,i) = 1;
                    end
                end
                obj.CosCal = z_v * z_d * z_v';
                A = obj.CosCal;
                obj.isCosCal = true;
            else
                A = obj.CosCal;
            end
        end
        % TODO: fix using multiple sin functions. should use only one with
        % Arnoldi it
        function A = sin_hTOmega(obj)
            if ~obj.isSinCal
                z_v = obj.O2_V;
                z_d = obj.O2_D;
                h = obj.hT(end);
                for i = 1:length(z_v)
                    if(z_d(i,i) ~= 0)
                        z_d(i,i) = sin(h*sqrt(z_d(i,i)));
                    else
                        z_d(i,i) = 0;
                    end
                end
                obj.SinCal = z_v * z_d * z_v';
                A = obj.SinCal;
                obj.isSinCal = true;
            else
                A = obj.SinCal;
            end
        end
        
        function A = osin_hTOmega(obj)
            if ~obj.isOsinCal
                z_v = obj.O2_V;
                z_d = obj.O2_D;
                h = obj.hT(end);
                for i = 1:length(z_v)
                    if(z_d(i,i) ~= 0)
                        z_d(i,i) = sqrt(z_d(i,i))*sin(h*sqrt(z_d(i,i)));
                    else
                        z_d(i,i) = 0;
                    end
                end
                obj.OsinCal = z_v * z_d * z_v';
                A = obj.OsinCal;
                obj.isOsinCal = true;
            else
                A = obj.OsinCal;
            end
        end
        function A = sinc_hTOmega(obj)
            if ~obj.isSincCal
                z_v = obj.O2_V;
                z_d = obj.O2_D;
                h = obj.hT(end);
                for i = 1:length(z_v)
                    if(z_d(i,i) ~= 0)
                        z_d(i,i) = sinc(h/pi*sqrt(z_d(i,i)));
                    else
                        z_d(i,i) = 1;
                    end
                end
                obj.SincCal = z_v * z_d * z_v';
                A = obj.SincCal;
                obj.isSincCal = true;
            else
                A = obj.SincCal;
            end
        end
        
        function A = sinc2_hTOmega(obj)
            if ~obj.isSinc2Cal
                z_v = obj.O2_V;
                z_d = obj.O2_D;
                h = obj.hT(end);
                for i = 1:length(z_v)
                    if(z_d(i,i) ~= 0)
                        z_d(i,i) = (sinc(h*sqrt(z_d(i,i))))^2;
                    else
                        z_d(i,i) = 1;
                    end
                end
                obj.Sinc2Cal = z_v * z_d * z_v';
                A = obj.Sinc2Cal;
                obj.isSinc2Cal = true;
            else
                A = obj.Sinc2Cal;
            end
        end
        
        function A = phi_hTOmega(obj)
            % for options see p.256 in M. Hochbruck and A. Ostermann
            % "exponential integrators"
            A = obj.sinc_hTOmega; % use phi = sinc
        end
        
        function A = psi_hTOmega(obj)
            % for options see p.256 in M. Hochbruck and A. Ostermann
            % "exponential integrators"
            A = obj.sinc2_hTOmega; % use psi = sinc^2
            
        end
        
        function A = psi0_hTOmega(obj)
            A = obj.sinc_hTOmega * obj.cos_hTOmega;
        end
        
        function A = psi1_hTOmega(obj)
            A = obj.sinc_hTOmega;
        end
        
        function s = g(obj,q)
            %  use the nonlinear part of M\f as g
            s = obj.M\obj.HD(obj.T,q) + obj.HOmega2(obj.T,q)*q;
        end
        function s = Mg(obj,q)
        %  use the nonlinear part of f as Mg
            s = obj.HD(obj.T,q) + obj.M*obj.HOmega2(obj.T,q)*q;
        end
        
        function O = HOmega2(obj,t,q) % function handle to Omega2
            O = -full(obj.M\obj.HDJ(t,q));
        end
    end
    
end
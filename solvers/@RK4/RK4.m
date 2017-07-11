classdef RK4 < FirstOrderIVP
    %RK4 Solving a first-order IVP by 4th order RK
    %   Detailed explanation goes here
    properties
    end

    methods
        function obj = RK4(Hf, IC)
            if nargin == 0
                error('Please provide the function handle to the derivative and the initial condition')
            end
            obj = obj@FirstOrderIVP(Hf, IC);

        end
        
        function sol = step(obj,h)
            % 
            assert(isa(obj.Hf,'function_handle'));
            obj.state = [obj.state; zeros(1,size(obj.state,2))];
            t = obj.T;
            u = obj.state(obj.CurrentIndex,:)';
            k1 = h*obj.Hf(t,u);
            k2 = h*obj.Hf(t+h/2, u+k1/2);
            k3 = h*obj.Hf(t+h/2, u+k2/2);
            k4 = h*obj.Hf(t+h, u+k3);
            obj.state(obj.CurrentIndex + 1,:) = u + k1/6 + k2/3 + k3/3 + k4/6;
            
            obj.CurrentIndex = obj.CurrentIndex + 1;
            obj.T = obj.T + h;
            obj.hT = [obj.hT h];
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
    
    methods(Static)
        test_RK4() % testing RK4 with some simple ode 
        test_RK4_sine_Gordon()
        test_RK4_nonlinear_wave()
    end
    
end


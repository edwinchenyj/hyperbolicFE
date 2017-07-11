classdef ForwardEuler < FirstOrderIVP
    %ForwardEuler Solving a first-order IVP by ForwardEuler
    %   Detailed explanation goes here
    properties
    end

    methods
        function obj = ForwardEuler(Hf, IC)
            if nargin == 0
                error('Please provide the function handle to the derivative and the initial condition')
            end
            obj = obj@FirstOrderIVP(Hf, IC);
%             assert(isa(obj.Hf,'function_handle'));
%             
%             obj.Hf = Hf;
%             obj.IC = IC;
%             obj.CurrentIndex = 1;
%             obj.state(1,:) = IC;
%             obj.hT = [];
%             obj.T = [];
        end
        
        function sol = step(obj,h)
            assert(isa(obj.Hf,'function_handle'));
            obj.state = [obj.state; zeros(1,size(obj.state,2))];
            u = obj.state(obj.CurrentIndex,:)'; % convert to a column vector
            obj.state(obj.CurrentIndex + 1,:) = (u + h * obj.Hf(obj.T,u))';
            
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
        test_forward_euler()
    end
end


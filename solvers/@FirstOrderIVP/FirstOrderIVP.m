classdef FirstOrderIVP < handle
    %FirstOrderIVP Describe a first order IVP and record the solution
    %states
    %   Detailed explanation goes here
    %   u' = f(u,t)
    
    properties (SetAccess = protected, GetAccess = protected)
        IC; % initial condition (row vector)
        hT; % a column vector storing all the time steps taken
        T = 0; % current time
        Hf; % the function handle to f(u,t)
        state; % the state of the IVP at each time step
        
        CurrentIndex = 1; % the index to the current time and state;
    end
    
    methods
        function obj = FirstOrderIVP(Hf, IC)
            if nargin == 0
                warning('Please provide the function handle to the derivative and the initial condition')
            else
                assert(isa(Hf,'function_handle'));
                
                obj.Hf = Hf;
                obj.IC = IC;
                obj.state(1,:) = IC;
                obj.hT = [];
                obj.T = 0; % default starting time at 0

            end
        end
        function CurrentIndex = GetCurrentIndex(obj)
            CurrentIndex = obj.CurrentIndex;
        end
        
        function CurrentState = GetCurrentState(obj)
            CurrentState = obj.state(obj.CurrentIndex,:);
        end
        
        function CurrentTime = GetCurrentTime(obj)
            CurrentTime = obj.hT(obj.CurrentIndex);
        end
        
        function State = GetState(obj,index)
            State = obj.state(index,:);
        end
        
        function StepSize = GetStepSize(obj,index)
            StepSize = obj.hT(index);
        end
        
        function ic = GetIC(obj)
            ic = obj.IC;
        end
        LoadStates(obj,u,t) % load the states previously computed
        
        function SaveStates(obj) % Save the entire ODE states into 'state' and 'hT'
            classname = class(obj);
            i = 1;
            filename = [classname num2str(i)];
            while exist([filename '.mat'], 'file') == 2
                i = i + 1;
                filename = [classname num2str(i)];
            end
%             firstOrderIVP = obj;
            save(filename);
            
        end
    end
    
    methods (Abstract = true)
        sol = solve(obj, T) % evolve the state along the time in T
        sol = step(obj,h) % advance by time step h and return the new state
    end
    

    
end


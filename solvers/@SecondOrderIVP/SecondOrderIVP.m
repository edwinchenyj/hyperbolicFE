classdef SecondOrderIVP < handle
    %SecondOrderIVP Describe a second order IVP and record the solution
    %states
    %   Detailed explanation goes here
    %   q'' + M\B(q,q') * q' = M\D(t,q) = -M\K(q) * q + M\f(t,q) 
    
    properties (SetAccess = protected, GetAccess = protected)

        hT; % a column vector storing all the time steps taken
        T = 0; % current time
        HD; % the function handle to D(t,q)
        HDJ; % the function handle to DJ(t,q), the Jacobain of D = -K + linearized(f)
        HB; % the function handle to B(q,q')
        state; % the full state of the IVP at each time step
        IC; % initial condition (row vector)
        M; % mass matrix
        CurrentIndex = 1; % the index to the current time and state;
    end
    
    methods
        function obj = SecondOrderIVP(HD, HDJ, IC, M, varargin)
            % varargin{1} should contain the function handle to B(u,u') =
            % HB
            if nargin == 0
                warning('Please provide the function handle to the derivative and the initial condition')
            else
                assert(isa(HD,'function_handle'));
                
                obj.HD = HD;
                obj.HDJ = HDJ;
                obj.IC = IC;
                obj.M = M;
                    obj.HB = varargin{1}
                obj.state(1,:) = IC;
                obj.hT = [];
                obj.T = 0; % default starting time at 0
                
                if nargin > 4;
                end
                
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
        
    end
    
    methods (Abstract = true)
        sol = solve(obj, T) % evolve the state along the time in T
        sol = step(obj,h) % advance by time step h and return the new state
    end
    

    
end


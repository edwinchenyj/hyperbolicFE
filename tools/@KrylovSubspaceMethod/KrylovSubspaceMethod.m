classdef (Abstract = true) KrylovSubspaceMethod < handle
    %KrylovSubspaceMethod A superclass of krylov subspace methods
    %implemented in a matrix-free fashioin
    %   Detailed explanation goes here
    
    properties (Access = protected)
        HMV; % the function handle to the matrix-vector product function
        
        tol = 1e-6;
    end
    
    methods
        function obj = KrylovSubspaceMethod(HMV)
            if nargin > 0
                assert(isa(HMV,'function_handle'));
                obj.HMV = HMV;
            end
        end
        
        function SetHMV(obj, HMV)
            obj.HMV = HMV;
        end        
        
        function SetWarmStart(obj, x0)
            obj.isWarmStart = true;
            obj.x0 = x0;
        end
    end
    
end


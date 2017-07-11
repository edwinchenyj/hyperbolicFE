classdef Newton < handle
    % A class of nonlinear solver using Newton's method
    % i.e. solves F(y) = 0
    properties
        HF; 
        % the function handle to F(y,varargin). Notice F may have other
        % input parameters
        HDF;
        % the Jacobian of the input function. 
        isMatrixFree = false;
        init;
        % initial guess of the solution
        N;
        HDFv;
        % the function handle to the Jacobian-vector product. Used when
        % using matrix-free calculation
        MaxIT = 400;
    end
    methods
        function obj = Newton(HF, init, varargin)
%           varargin list if matrix free:
%           {HDF, true}
%           varargin list if not matrix free:
%           {HDF, false}
            obj.HF = HF;
            obj.init = init;
            obj.N = length(init);
            if nargin > 2
                obj.HDF = varargin{1};
                if varargin{2}
                    obj.isMatrixFree = true;
                end
            end
                
        end
        
        function [flag, sol] = solve(obj)
            D = Inf(obj.N,1);
            sol = obj.init;
            it = 0;
            if obj.isMatrixFree
                % TODO: support matrix free
            else
                % inexact Newton: fixing the Jacobian and the nonlinear
                % part for each Newton iterations
                % J = obj.HJ(t, OldState);
                while (norm(D) > 1e-8) && (it < obj.MaxIT)
                    % for full newton, update J in each iteration
                    J = obj.HDF(sol);
                    
                    D = -J\obj.HF(sol);
                    sol = sol + D;
                    it = it + 1;
                end
            end
            if it == obj.MaxIT
                % non-converging flag
                flag = true;
            else
                flag = false;
            end
        end
    end
    
    methods(Static)
        test_newton;
    end
end
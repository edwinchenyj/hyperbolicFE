classdef CG < KrylovSubspaceMethod
    %CG Solve A\b using conjugate gradient in a matrix-free fashion. Must
    %provide the function handle to the function of the matrix-vector
    %product A*v
    % TODO: add assertion for SPD
    %   Detailed explanation goes here
    
    properties
        tol = 1e-6;
        MaxIT = 40;
        isWarmStart = false;
        x0; % initial guess
    end
    
    methods
        function obj = CG(varargin)
            % have the option of warm-starting the solver by passing in an
            % optional initial guess x0
            if nargin > 0
                HMV = varargin{1};
                obj.HMV = HMV;
            if nargin == 2
                obj.isWarmStart = true;
                obj.x0 = varargin{2};
            end
            end
        end
        
        
        function [sol, flag] = solve(obj,b)
            if obj.isWarmStart
                xkm1 = obj.x0; % xk-1
                xk = obj.x0;
            else
                xkm1 = zeros(size(b));
                xk = zeros(size(b));
            end
            rk = b - obj.HMV(xkm1); % r0 = b - A * x0
            deltakm1 = rk'*rk; 
            deltak = deltakm1;
            bd = b' * b; % b_d = < b, b >
            pk = rk; % p0 = r0
            
            k = 0;
            flag = true;
            while (deltak > obj.tol^2 * bd)
                if k > obj.MaxIT
                    flag = false;
                    break;
                end
                sk = obj.HMV(pk);
                alphak = deltak / (pk' * sk);
                xk = xkm1 + alphak * pk; % or xk+1 = xk + alphak * pk
                rk = rk - alphak * sk;
                deltak = rk' * rk;
                pk = rk + deltak/deltakm1 * pk;
                deltakm1 = deltak;
                xkm1 = xk;
                k = k + 1;
            end
            sol = xk;
        end
    end
    
    methods(Static)
        test_cg()
    end
    
end


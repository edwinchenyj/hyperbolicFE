classdef GMRES < KrylovSubspaceMethod
    %GMRES Solve A\b using GMRES with modified Gram-Schmidt in a matrix-free fashion. Must
    %provide the function handle to the function of the matrix-vector
    %product A*v
    % TODO: add termination condtion for GMRES, see Yousef Saad section
    % 6.5.3 at page 167
    %   Detailed explanation goes here
    
    properties
        MaxD = Inf;
        MaxIT = 40;
        isWarmStart = false;
        x0; % initial guess
    end
    
    methods
        function obj = GMRES(varargin)
            % have the option of warm-starting the solver by passing in an
            % optional initial guess x0
            if nargin > 0
                HMV = varargin{1};
                obj.HMV = HMV;
            if nargin == 2
                obj.isWarmStart = true;
                obj.x0 = varargin{2};
                obj.MaxIT = min([obj.MaxIT length(obj.x0)]);
            end
            end
        end
                
        function [sol] = solve(obj,b)
            
            if obj.isWarmStart
                x_init = obj.x0; % xk-1
            else
                x_init = zeros(size(b));
            end
            r0 = b - obj.HMV(x_init); % r0 = b - A * x0
            beta = sqrt(r0'*r0); 
            v1 = r0 / beta;
            
            m = min(obj.MaxD, length(b));
            h = zeros(m+1);
            v = zeros(length(v1),m); v(:,1) = v1;
            w = zeros(size(v));
            for j = 1:m
                w(:,j) = obj.HMV(v(:,j));
                for i = 1:j
                    h(i,j) = w(:,j)' * v(:,i);
                    w(:,j) = w(:,j) - h(i,j)*v(:,i);
                end
                h(j+1,j) = sqrt(w(:,j)' * w(:,j));
                if abs(h(j+1,j)) < 1e-3
                    m = j;
                    break;
                end
                v(:,j+1) = w(:,j)/h(j+1,j);
            end
            H = h(1:m+1,1:m);
            e1 = zeros(m+1,1);
            e1(1) = 1;
            y = H\(beta* e1);
            sol = x_init + v(:,1:m) *y;
        end
    end
    
    methods(Static)
        test_gmres()
    end
    
end


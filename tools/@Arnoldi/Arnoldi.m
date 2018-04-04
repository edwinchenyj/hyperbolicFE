classdef Arnoldi < KrylovSubspaceMethod
    % Arnoldi's method is an orthogonal projection on to the krylov
    % subspace
    % Input: 
    % 1. The function handle to the function of the matrix-vector
    % product A*v
    % 2. (varargin) The input v
    % 3. (varargin) the dimension of the krylov subspace m
    %   Detailed explanation goes here
    
    properties
        v; % initial vector
        m = 10; % dimension of the krylov subspace
    end
    
    methods
        function obj = Arnoldi(varargin)
            % TODO: formalize the input structure
            if nargin > 0
                obj.HMV = varargin{1};
                if nargin > 1
                    obj.v = varargin{2};
                    if nargin > 2
                        obj.m = varargin{3};
                    end
                end
            end
        end
        
        
        function [Vmp1, HmBar, D] = ConstructBasis(obj)
            % Construct the basis using MGS
            % Output: Vm, HmBar: projection
            %       D = dimension
            % TODO: add support for Householder (for better stability)
            
            Vmp1 = zeros(length(obj.v),obj.m);
            Vmp1(:,1) = normc(real(obj.v(:)));
            w = zeros(length(obj.v),obj.m);
            HmBar = zeros(obj.m+1,obj.m);
            D = 1;
            for j = 1:obj.m
                w(:,j) = obj.HMV(Vmp1(:,j));
                for i = 1:j
                    HmBar(i,j) = w(:,j)' * Vmp1(:,i);
                    w(:,j) = w(:,j) - HmBar(i,j) * Vmp1(:,i);
                end
                HmBar(j+1, j) = norm(w(:,j));
                if abs(HmBar(j+1, j)) < 1e-6
                    
                    break;
                end
                Vmp1(:,j+1) = w(:,j)/ HmBar(j+1,j);
                D = D + 1;
            end
        end
    end

    methods(Static)
        test_arnoldi()
    end
    
end


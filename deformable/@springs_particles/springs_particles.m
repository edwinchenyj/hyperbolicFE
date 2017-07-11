classdef springs_particles < handle
    % A system of springs connected particles
    %   Detailed explanation goes here
    
    properties
        n; % number of particles
        m; % number of springs
        k; % stiffness
        k_damping; % damping constants
        q; % position state
        v; % velocity state
        springs;
        % each row contains the indices of the particles connected by that
        % spring
        l; % rest length of the springs
        mass; % mass
        mass_matrix;
        k_matrix; % the stiffness matrix
    end
    
    methods
        function obj = springs_particles(q,v,springs,k,k_damping,l,mass)
            obj.n = length(q)/3;
            obj.m = length(l);
            obj.q = q;
            obj.v = v;
            obj.springs = springs;
            obj.k = k;
            obj.k_damping = k_damping;
            obj.l = l;
            obj.mass = mass;
            obj.mass_matrix = diag(mass);
            obj.k_matrix = sparse(3*obj.n,3*obj.n);
        end
        
        function f = force(obj)
            
            
            f = zeros(3 * obj.n , 1);
            % f(z_ind) = -9.8 * obj.mass(z_ind);
            for sp = 1:obj.m
                i = obj.springs(sp,1);
                j = obj.springs(sp,2);
                
                vector_l = obj.q(3*(i-1) + 1 : 3*i) - obj.q(3*(j-1) + 1 : 3*j);
                norm_l = norm(vector_l);
                % strain of the spring
                strain = vector_l / norm_l * (norm_l - obj.l(sp));
                % spring forces
                fs = -obj.k(sp) *  strain;
                % damping forces
%                 delta_v = obj.v(3*(i-1) + 1 : 3*i) - obj.v(3*(j-1) + 1 : 3*j);
%                 fd = obj.k_damping(sp) *  dot(delta_v, vector_l)    / norm_l;
                fd = zeros(size(fs));

                f(3*(i-1) + 1 : 3*i) = f(3*(i-1) + 1 : 3*i) +fs +fd;
                f(3*(j-1) + 1 : 3*j) = f(3*(j-1) + 1 : 3*j) -fs -fd;
            end
        end
        
        function update(obj,q,v)
            obj.q = q;
            obj.v = v;
        end
        
        function K = stiffness_matrix(obj)
            % the stiffness matrix
            n = obj.n;
            q = obj.q;
            v = obj.v;
            k = obj.k;
            l = obj.l;
            I = eye(3);
            obj.k_matrix = sparse(3*n,3*n);
            
            for sp = 1:obj.m
                ks = k(sp);
                l0 = l(sp);
                i = obj.springs(sp,1);
                j = obj.springs(sp,2);
                vector_l = q(3*(i-1) + 1 : 3*i) - q(3*(j-1) + 1 : 3*j);
                norm_l = norm(vector_l);
                k_ii = -ks * (I - l0/norm_l * ( I - vector_l * vector_l'/(norm_l^2)));
                k_jj = k_ii;
                k_ij = - k_ii;
                k_ji = k_ij;
                
                obj.k_matrix(3*(i-1) + 1 : 3*i, 3*(i-1) + 1 : 3*i) = obj.k_matrix(3*(i-1) + 1 : 3*i, 3*(i-1) + 1 : 3*i) + k_ii;
                obj.k_matrix(3*(j-1) + 1 : 3*j, 3*(i-1) + 1 : 3*i) = obj.k_matrix(3*(j-1) + 1 : 3*j, 3*(i-1) + 1 : 3*i) + k_ji;
                obj.k_matrix(3*(i-1) + 1 : 3*i, 3*(j-1) + 1 : 3*j) = obj.k_matrix(3*(i-1) + 1 : 3*i, 3*(j-1) + 1 : 3*j) + k_ij;
                obj.k_matrix(3*(j-1) + 1 : 3*j, 3*(j-1) + 1 : 3*j) = obj.k_matrix(3*(j-1) + 1 : 3*j, 3*(j-1) + 1 : 3*j) + k_jj;
            
            end
            
            K = obj.k_matrix;
            
        end
        
        
    end
    
    methods (Static)
        test_directional_derivative()
    end
    
end


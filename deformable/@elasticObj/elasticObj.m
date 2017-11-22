classdef (Abstract) elasticObj < handle

    % elasticTriObj An elastic object represented by a triangle mesh
    %   Detailed explanation goes here
    
    
    properties
        is_polyfit = false;
        is_DAC = false;
        N; % #nodes
        NT; % #elements
        
        node; % nodal position in world (N by 2)
        nodeM; % nodal position in material (N by 2)
        elem; % element label (NT by 3): [N1 N2 N3]
        
        x; % position vector in world (2N by 1)
        v; % velocity vector in world (2N by 1)
        X; % position vector in material, AKA rest pose (2N by 1)
        
        Ds; % some type of discrete gradient per element (2NT by 2) storing 2 edges of the elements in world
        Dm; % some type of discrete gradient per element (2NT by 2) storing 2 edges of the elements in material
        DmINV; % inverse of Dm (2NT by 2)
        W; % undeformed volume of each element (NT by 1)
        T; % mapping vectorized nodal position in a tri to its vectorized deformation gradient (4NT by 6)
        % definition: vec(F) = T * vec(x), or vec(dF) = T * vec(dx)
        M % mass matrix
        K0 = [];
        
        ii; % sparse structure
        jj; % sparse structrue for the stiffness matrix
        
        F; % deformation gradient (2NT by 2) from F*Dm = Ds = xG and F = xG(Dm)^(-1)
        FINV % inverse of F (2NT by 2)
        
        U; % left singular vectors for each tri
        V; % right singular vectors for each tri
        S; % principle stretches for each tri (singular values)
        R; % rotation for each tri
        
        % TODO: may want to add J = det(F), or I1, I2, and I3, the
        % invariances if the deformation gradient tensor
        C; % green strain
        normC;
        
        f; % elastic force under the current deformation
        
        isCorotatedLinear = false;
        material_type;
        
        mu; % An array of mu to the elements 
        lambda; % An array of lambda to the elements 
        rho; % An array of rho to the elements 
        
        % TODO: need to get rid of this
        Y;
        P;
        Rho;
        
        % parameters for Rayleigh damping
        a = 0;
        b = 0;
        
        gravity_on = false;
        externalGravity;
        
        % for ghost sub-points
        sub_objs;
        
        % handle to the axis for visualization
        vis_handle;
        
        % force approximation for polyfit
        polyfit_force_approx;
        % polynomial
        polyfit_p;
        % if the simulating start from rest state, not extra computation
        % needed
        from_rest = true;
    end
    
    methods (Abstract)
        
        ElasticForce(obj)
        
        simple_vis(obj)
        
        ElasticForceDifferential(obj)
        
        StiffnessMatrix(obj)
        
        SetMaterial(obj)
        
        totalEnergy(obj)
        
%         calculateGravity(obj)
        
        mesh_quality(obj)
        
        Stress(obj)
        
        StressDifferential(obj)
    end
    
    methods

        function ha = init_vis(obj,varargin)
            switch nargin
                case 1
                    ha = axes;
                    obj.vis_handle = ha;
                    hold(obj.vis_handle,'on');
                case 2 % case where the axes is provided
                    obj.vis_handle = varargin{1};
                    hold(obj.vis_handle,'on');
            end
        end
        function K = restStiffness(obj)
            assert(isequal(obj.x,obj.X))
            obj.K0 = obj.StiffnessMatrix;
            K = obj.K0;
        end
        
        function calculateGravity(obj)
            % by default the is no gravity. if need to simulate gravity,
            % need to call this function explicitly before hand
            if obj.gravity_on
                obj.externalGravity = repmat(obj.gravity,obj.N,1);
            else
                obj.externalGravity = zeros(size(obj.X,1),1);
            end
        end



    end

end
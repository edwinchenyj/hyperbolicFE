classdef elasticObj < handle
    %elasticObj An elastic object represented by a triangular mesh
    %   Detailed explanation goes here
    
    properties
        
        gravity = [0 0 0]';
        
        % to be removed
        C;
        dP;
        dFdx;
        dF;
    end
    
    properties (Constant)
        G = [1 0; 0 1; -1 -1]; % element-wise operator that maps element nodal position in to edge matrix Ds or Dm
        I2 = eye(2);
        I4 = eye(4);
        % commutation matrix
        K44 =   ...
            [1     0     0     0;...
            0     0     1     0;...
            0     1     0     0;...
            0     0     0     1];
    end
    properties (Access = private)
        N; % #nodes
        NT; % # triangular elements
        
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
        T; % mapping vectorized nodal position in a tet to its vectorized deformation gradient (4NT by 6)
        % definition: vec(F) = T * vec(x), or vec(dF) = T * vec(dx)
        
        F; % deformation gradient (2NT by 2) from F*Dm = Ds = xG and F = xG(Dm)^(-1)
        FINV % inverse of F (2NT by 2)
        % TODO: may want to add J = det(F), or I1, I2, the
        % invariances if the deformation gradient tensor
        
        f; % elastic force under the current deformation
        
        materialBlock = {}; % A cell array of indice to the elements of the same material type
        materialType = {}; % A cell array of type of material of each block
        materialBlockCount = 0; % number of different blocks of material
        elemMaterialType; % an array containing the material type of each element (NT by 1), 1 = neo-hookean, 2 = linear
        
        mu; % An array of mu to the elements 
        lambda; % An array of lambda to the elements 
        rho; % An array of rho to the elements 
        
        finalized = false; % A flag to check if the object is ready for simulation
        

    end
    
    methods
        function obj = elasticTetObj(varargin)
            % Constructor of elasticTetObj
            %   Take in an undeformed tet mesh with
            %       nodeM: (undeformed) nodal position in material (#nodes
            %       by 2)
            %       elem: element label (#element by 4) [N1 N2 N3]
            if nargin == 1 && isa(varargin{1},'elasticObj')
                % a deep-copy copy constructor
                oldObj = varargin{1};
                mc = metaclass(oldObj);
                metaPropertiesList = mc.PropertyList; 
                % get the meta.property from the meta.class from the
                % metaclass function, otherwise we may miss the private or
                % hidden properties
                for i = 1:length(metaPropertiesList)
                    if ~metaPropertiesList(i).Constant
                        % make sure we are not writting a Constant property
                        obj.(metaPropertiesList(i).Name) = oldObj.(metaPropertiesList(i).Name);
                    end
                end
            elseif nargin > 1
                input_nodeM = varargin{1};
                input_elem = varargin{2};
                assert(size(input_elem,2)== 3);
                assert(size(input_nodeM,2) == 2);
                
                obj.N = size(input_nodeM,1);
                obj.NT = size(input_elem,1);
                obj.elemMaterialType = zeros(obj.NT,1);
                obj.mu = zeros(obj.NT,1);
                obj.lambda = zeros(obj.NT,1);
                obj.rho = zeros(obj.NT,1);
                
                obj.nodeM = input_nodeM;
                obj.elem = input_elem;
                
                obj.X = reshape(input_nodeM',2*obj.N,1);
                obj.Dm = zeros(2*obj.NT,2);
                obj.DmINV = zeros(2*obj.NT,2);
                obj.W = zeros(obj.NT,1);
                obj.T = zeros(4*obj.NT, 12);
                for i = 1:obj.NT
                    T_node = input_nodeM(input_elem(i,:),:); % element nodal position in material space  (#nodes per elements by 2)
                    obj.Dm(2*(i-1)+1:2*i,:) = T_node'*obj.G;
                    obj.DmINV(2*(i-1)+1:2*i,:) = inv(T_node'*obj.G);
                    obj.W(i) =  -det(T_node'*obj.G)/6; % undeformed volume from matrix determinant
                    % using the negative determinant because 
                    obj.T(9*(i-1)+1:9*i,:) = kron((obj.G * obj.DmINV(2*(i-1)+1:2*i,:))', obj.I2); % definition: vec(F) = T * vec(x), or vec(dF) = T * vec(dx)
                end
            end
        end
        
        SetCurrentState(obj, x, v)
        
        SetMaterial(obj, Y, P, Rho, elem, type)
        
        function N = GetNNodes(obj)
            N = obj.N;
        end
        
        function state = GetX(obj)
            state = obj.x;
        end
        
        function v = GetV(obj)
            v = obj.v;
        end
        
        f = ElasticForce(obj)

        df = ElasticForceDifferential(obj, dx)
        
        K = StiffnessMatrix(obj)
        
        function finalized = isFinalized(obj)
            % return true is the object is ready for simulation
            finalized = obj.finalized;
        end
        
        function finalize(obj)
            % mark finalize true if all the required fields are filled
            % TODO need to check the current state and the undeformed state exist
            for t = 1:obj.NT
                assert(obj.elemMaterialType(t) ~= 0); % check all elements have been initialized
            end
            obj.finalized = true;
        end
    end
    
    methods (Access = private)
        P = Stress(obj, t)
        
        dP = StressDifferential(obj, t, dF)

    end
    methods (Static)
        test_matrix_free()
        test_stress()
        test_elastic_force()
        test_object_deep_copy()
        test_tet_orientation()
    end
end

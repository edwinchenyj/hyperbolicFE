classdef elasticTetObj < handle
    %elasticTetObj An elastic object represented by a tet mesh
    %   Detailed explanation goes here
    
    properties
        
        gravity = [0 0 9.81]';
        
        % to be removed
        
        dP;
        dFdx;
        dF;
    end
    
    properties (Constant)
        G = [1 0 0; 0 1 0; 0 0 1; -1 -1 -1]; % element-wise operator that maps element nodal position in to edge matrix Ds or Dm
        I3 = eye(3);
        I9 = eye(9);
        % commutation matrix
        K99 =   ...
            [1     0     0     0     0     0     0     0     0;...
            0     0     0     1     0     0     0     0     0;...
            0     0     0     0     0     0     1     0     0;...
            0     1     0     0     0     0     0     0     0;...
            0     0     0     0     1     0     0     0     0;...
            0     0     0     0     0     0     0     1     0;...
            0     0     1     0     0     0     0     0     0;...
            0     0     0     0     0     1     0     0     0;...
            0     0     0     0     0     0     0     0     1];
        
        IndK = [[1:3, [1:3] + 12, [1:3] + 24]; ...
            [1:3, [1:3] + 12, [1:3] + 24] + 36; ...
            [1:3, [1:3] + 12, [1:3] + 24] + 72; ...
            [1:3, [1:3] + 12, [1:3] + 24] + 108; ...
            ...
            [1:3, [1:3] + 12, [1:3] + 24] + 3; ...
            [1:3, [1:3] + 12, [1:3] + 24] + 36 + 3; ...
            [1:3, [1:3] + 12, [1:3] + 24] + 72 + 3; ...
            [1:3, [1:3] + 12, [1:3] + 24] + 108 + 3; ...
            ...
            [1:3, [1:3] + 12, [1:3] + 24] + 6; ...
            [1:3, [1:3] + 12, [1:3] + 24] + 36 + 6; ...
            [1:3, [1:3] + 12, [1:3] + 24] + 72 + 6; ...
            [1:3, [1:3] + 12, [1:3] + 24] + 108 + 6; ...
            ...
            [1:3, [1:3] + 12, [1:3] + 24] + 9; ...
            [1:3, [1:3] + 12, [1:3] + 24] + 36 + 9; ...
            [1:3, [1:3] + 12, [1:3] + 24] + 72 + 9; ...
            [1:3, [1:3] + 12, [1:3] + 24] + 108 + 9; ...
            ]; % fast indexing for the stiffness matrix (for optimizating the speed)
    end
    properties
        N; % #nodes
        NT; % #elements
        
        node; % nodal position in world (N by 3)
        nodeM; % nodal position in material (N by 3)
        elem; % element label (NT by 4): [N1 N2 N3 N4]
        
        x; % position vector in world (3N by 1)
        v; % velocity vector in world (3N by 1)
        X; % position vector in material, AKA rest pose (3N by 1)
        
        Ds; % some type of discrete gradient per element (3NT by 3) storing 3 edges of the elements in world
        Dm; % some type of discrete gradient per element (3NT by 3) storing 3 edges of the elements in material
        DmINV; % inverse of Dm (3NT by 3)
        W; % undeformed volume of each element (NT by 1)
        T; % mapping vectorized nodal position in a tet to its vectorized deformation gradient (9NT by 12)
        % definition: vec(F) = T * vec(x), or vec(dF) = T * vec(dx)
        M % mass matrix
        K0;
        
        ii; % sparse structure
        jj; % sparse structrue for the stiffness matrix
        
        F; % deformation gradient (3NT by 3) from F*Dm = Ds = xG and F = xG(Dm)^(-1)
        FINV % inverse of F (3NT by 3)
        
        U; % left singular vectors for each tet
        V; % right singular vectors for each tet
        S; % principle stretches for each tet (singular values)
        R; % rotation for each tet
        
        % TODO: may want to add J = det(F), or I1, I2, and I3, the
        % invariances if the deformation gradient tensor
        C; % green strain
        normC;
        
        f; % elastic force under the current deformation
        
        isCorotatedLinear = false;
        materialBlock = {}; % A cell array of indice to the elements of the same material type
        materialType = {}; % A cell array of type of material of each block
        materialBlockCount = 0; % number of different blocks of material
        elemMaterialType; % an array containing the material type of each element (NT by 1), 1 = neo-hookean, 2 = linear
        
        mu; % An array of mu to the elements 
        lambda; % An array of lambda to the elements 
        rho; % An array of rho to the elements 
        
        % TODO: need to get rid of this
        Y;
        Rho;
        

    end
    
    methods
        function obj = elasticTetObj(varargin)
            % Constructor of elasticTetObj
            %   Take in an undeformed tet mesh with
            %       nodeM: (undeformed) nodal position in material (#nodes
            %       by 3)
            %       elem: element label (#element by 4) [N1 N2 N3 N4]
            if nargin == 1 && isa(varargin{1},'elasticTetObj')
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
                assert(size(input_elem,2)== 4);
                assert(size(input_nodeM,2) == 3);
                
                obj.N = size(input_nodeM,1);
                obj.NT = size(input_elem,1);
                obj.elemMaterialType = zeros(obj.NT,1);
                obj.mu = zeros(obj.NT,1);
                obj.lambda = zeros(obj.NT,1);
                obj.rho = zeros(obj.NT,1);
                
                obj.nodeM = input_nodeM;
                obj.elem = input_elem;
                
                obj.X = reshape(input_nodeM',3*obj.N,1);
                obj.Dm = zeros(3*obj.NT,3);
                obj.DmINV = zeros(3*obj.NT,3);
                obj.F = zeros(3*obj.NT,3);
                obj.U = zeros(3*obj.NT,3);
                obj.V = zeros(3*obj.NT,3);
                obj.S = ones(3*obj.NT,1);
                obj.R = zeros(3*obj.NT,3);
                obj.C = zeros(3*obj.NT,3);
                obj.normC = zeros(obj.NT,1);
                obj.M = sparse(size(obj.X,1),size(obj.X,1));
                obj.W = zeros(obj.NT,1);
                obj.T = zeros(9*obj.NT, 12);
                index = 1;
                for i = 1:obj.NT
                    T_node = input_nodeM(input_elem(i,:),:); % element nodal position in material space  (#nodes per elements by 3)
                    obj.Dm(3*(i-1)+1:3*i,:) = T_node'*obj.G;
                    obj.DmINV(3*(i-1)+1:3*i,:) = inv(T_node'*obj.G);
                    obj.W(i) =  -det(T_node'*obj.G)/6; % undeformed volume from matrix determinant
                    % using the negative determinant because of orientation
                    obj.T(9*(i-1)+1:9*i,:) = kron((obj.G * obj.DmINV(3*(i-1)+1:3*i,:))', obj.I3); % definition: vec(F) = T * vec(x), or vec(dF) = T * vec(dx)
                    for e = input_elem(i,:)
                        for mi = (e-1)*3+1:e*3
                            obj.M(mi,mi) = obj.M(mi,mi)+obj.W(i)/4;
                        end
                    end
                    
                    for ti = 1:4
                        for tj = 1:4
                            obj.ii(index:index+8) = repmat((3*(obj.elem(i,ti)-1)+1:3*obj.elem(i,ti))',3,1);
                            obj.jj(index:index+8) = reshape(repmat(3*(obj.elem(i,tj)-1)+1:3*obj.elem(i,tj),3,1),9,1);
                            index = index + 9;
                        end
                    end
                    
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
        
        function K = restStiffness(obj)
            assert(isequal(obj.x,obj.X))
            obj.K0 = obj.StiffnessMatrix;
            K = obj.K0;
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
        
        energy = totalEnergy(obj)
        
        function C = greenStrain(obj)
            for i = 1:obj.NT
                    obj.C(3*(i-1)+1:3*i,:) = 1/2 * (obj.F(3*(i-1)+1:3*i,:)'*obj.F(3*(i-1)+1:3*i,:)-eye(3));
            end
            C = obj.C;
        end
        
        function normC = greenStrainNorm(obj)
            obj.greenStrain;
            for i = 1:obj.NT
                    obj.normC(i) = norm(obj.C(3*(i-1)+1:3*i,:),'fro');
            end
            normC = obj.normC;
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
        test_elastic_force_2()
        test_object_deep_copy()
        test_tet_orientation()
    end
end

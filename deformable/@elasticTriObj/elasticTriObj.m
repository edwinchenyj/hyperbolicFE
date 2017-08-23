classdef elasticTriObj < handle

    % elasticTriObj An elastic object represented by a triangle mesh
    %   Detailed explanation goes here
    
    properties
        
        % to be removed
        Dim = 2;
        dP;
        dFdx;
        dF;
    end
    
    properties (Constant)

        G = [1 0; 0 1; -1 -1];
        I2 = eye(2);
        I4 = eye(4);
        % commutation matrix
        K44 = [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1];
        
        
        IndK = [[1:2, [1:2] + 6]; ...
            [1:2, [1:2] + 6] + 12; ...
            [1:2, [1:2] + 6] + 24; ...
            [1:2, [1:2] + 6] + 36; ...
            ...
            [1:2, [1:2] + 6] + 2; ...
            [1:2, [1:2] + 6] + 12 + 2; ...
            [1:2, [1:2] + 6] + 24 + 2; ...
            [1:2, [1:2] + 6] + 36 + 2; ...
            ...
            [1:2, [1:2] + 6] + 4; ...
            [1:2, [1:2] + 6] + 12 + 4; ...
            [1:2, [1:2] + 6] + 24 + 4; ...
            [1:2, [1:2] + 6] + 36 + 4; ...
            ]; % fast indexing for the stiffness matrix (for optimizating the speed)
    end
    properties
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
        K0;
        
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
        function obj = elasticTriObj(varargin)
            % Constructor 
            %   Take in an undeformed tri mesh with
            %       nodeM: (undeformed) nodal position in material (#nodes
            %       by 3)
            %       elem: element label (#element by 3) [N1 N2 N3]
            if nargin == 1 && isa(varargin{1},'elasticTriObj')
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
                obj.F = zeros(2*obj.NT,2);
                obj.U = zeros(2*obj.NT,2);
                obj.V = zeros(2*obj.NT,2);
                obj.S = ones(2*obj.NT,1);
                obj.R = zeros(2*obj.NT,2);
                obj.C = zeros(2*obj.NT,2);
                obj.normC = zeros(obj.NT,1);
                obj.M = sparse(size(obj.X,1),size(obj.X,1));
                obj.W = zeros(obj.NT,1);
                obj.T = zeros(4*obj.NT, 6);
                index = 1;
                for i = 1:obj.NT
                    T_node = input_nodeM(input_elem(i,:),:); % element nodal position in material space  (#nodes per elements by 3)
                    obj.Dm(2*(i-1)+1:2*i,:) = T_node'*obj.G;
                    obj.DmINV(2*(i-1)+1:2*i,:) = inv(T_node'*obj.G);
                    obj.W(i) =  -det(T_node'*obj.G)/2; % undeformed volume from matrix determinant
                    % using the negative determinant because of orientation
                    obj.T(4*(i-1)+1:4*i,:) = kron((obj.G * obj.DmINV(2*(i-1)+1:2*i,:))', obj.I2); % definition: vec(F) = T * vec(x), or vec(dF) = T * vec(dx)
                    
%                     % simple mass lumping by distributing equally
%                     % probably not the most accurate one...
%                     for e = input_elem(i,:)
%                         for mi = (e-1)*2+1:e*2
%                             obj.M(mi,mi) = obj.M(mi,mi)+obj.W(i)/3;
%                         end
%                     end
%                     
                    for ti = 1:3
                        for tj = 1:3
                            obj.ii(index:index+3) = repmat((2*(obj.elem(i,ti)-1)+1:2*obj.elem(i,ti))',2,1);
                            obj.jj(index:index+3) = reshape(repmat(2*(obj.elem(i,tj)-1)+1:2*obj.elem(i,tj),2,1),4,1);
                            index = index + 4;
                        end
                    end
                    
                end
            end
            
        end
        
        SetCurrentState(obj, x)
        
        SetMaterial(obj, Y, P, Rho, elem, type)
        
        function N = GetNNodes(obj)
            N = obj.N;
        end
        
        function state = GetX(obj)
            state = obj.x;
        end
        
        function K = restStiffness(obj)
            assert(isequal(obj.x,obj.X))
            obj.K0 = obj.StiffnessMatrix;
            K = obj.K0;
        end
        
        f = ElasticForce(obj)

        df = ElasticForceDifferential(obj, dx)
        
        K = StiffnessMatrix(obj)
        
        energy = totalEnergy(obj)
        
        function C = greenStrain(obj)
            for i = 1:obj.NT
                    obj.C(2*(i-1)+1:2*i,:) = 1/2 * (obj.F(2*(i-1)+1:2*i,:)'*obj.F(2*(i-1)+1:2*i,:)-eye(2));
            end
            C = obj.C;
        end
        
        function normC = greenStrainNorm(obj)
            obj.greenStrain;
            for i = 1:obj.NT
                    obj.normC(i) = norm(obj.C(2*(i-1)+1:2*i,:),'fro');
            end
            normC = obj.normC;
        end
        
        mesh_quality(obj)
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
        test_tri_orientation()
        test_mass()
        test_modes()
        test_modes_free()
        test_rect_free()
        test_bending()
    end

end
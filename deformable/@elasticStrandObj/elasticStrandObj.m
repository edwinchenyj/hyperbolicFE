classdef elasticStrandObj < handle
    %elasticObj An elastic object represented by a triangular mesh
    %   Detailed explanation goes here
    
    properties
        
        
        % to be removed
        C;
        dP;
        dFdx;
        dF;
    end
    
    properties (Constant)
        G = [-1; 1]; % element-wise operator that maps element nodal position in to edge matrix Ds or Dm

    end
    properties 
        N; % # of nodes
        NS; % # of strand elements
        
        node; % nodal position in world (N by 3)
        nodeM; % nodal position in material space (N by 3)
        elem; % element label (NS by 3): [N1 N2 N3]
        
        x; % position vector in world (3N by 1)
        v; % velocity vector in world (3N by 1)
        X; % position vector in material, AKA rest pose (3N by 1). For simple 1D strands, this will only be used to calculate the rest length
        
        Ds; % some type of discrete gradient per element (3NS by 1) storing 2 edges of the elements in world
        Dm; % some type of discrete gradient per element (3NS by 1) storing 2 edges of the elements in material
        W; % undeformed length of each element (NS by 1)
        T; % mapping vectorized nodal position in a strand to its vectorized deformation gradient (3NS by 6)
        % definition: vec(F) = T * vec(x), or vec(dF) = T * vec(dx)
        
        F; % deformation gradient (NS by 1) from F*Dm = Ds = xG and F = xG(Dm)^(-1)
        FINV % inverse of F (NS by 1)
        % TODO: may want to add J = det(F), or I1, I2, the
        % invariances if the deformation gradient tensor
        
        f; % elastic force under the current deformation
        
        materialBlock = {}; % A cell array of indice to the elements of the same material type
        materialType = {}; % A cell array of type of material of each block
        materialBlockCount = 0; % number of different blocks of material
        elemMaterialType; % an array containing the material type of each element (NS by 1), 1 = neo-hookean, 2 = linear
        
        ks; % An array of ks to the elements 
        rho; % An array of rho to the elements 
        M; % mass matrix
        K0; % rest stiffness
        
        ii; % 'sparse structure'
        jj; % 'sparse structure'
        
        finalized = false; % A flag to check if the object is ready for simulation
        
        

    end
    
    methods
        function obj = elasticStrandObj(varargin)
            % Constructor of elasticTetObj
            %   Take in an undeformed tet mesh with
            %       nodeM: (undeformed) nodal position in material (#nodes
            %       by 3)
            %       elem: element label (#element by 2) [N1 N2]
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
                assert(size(input_elem,2)== 2);
                assert(size(input_nodeM,2) == 3);
                
                obj.N = size(input_nodeM,1);
                obj.NS = size(input_elem,1);
                obj.elemMaterialType = zeros(obj.NS,1);
                obj.ks = zeros(obj.NS,1);
                obj.rho = zeros(obj.NS,1);
                
                obj.nodeM = input_nodeM;
                obj.elem = input_elem;
                
                obj.X = reshape(input_nodeM',3*obj.N,1);
                obj.Dm = zeros(3*obj.NS,1);
                obj.W = zeros(obj.NS,1);
                obj.T = zeros(obj.NS, 6);
                obj.M = sparse(3*obj.N, 3*obj.N);
                index = 1;
                for i = 1:obj.NS
                    S_node = input_nodeM(input_elem(i,:),:); % element nodal position in material space  (#nodes per elements by 3)
                    obj.Dm(3*(i-1)+1:3*i,:) = S_node'*obj.G;
                    obj.W(i) =  norm(S_node'*obj.G); % undeformed length from matrix determinant
                    obj.T(3*(i-1)+1:3*i,:) = [-eye(3) eye(3)]; % definition: vec(F) = T * vec(x), or vec(dF) = T * vec(dx)
                    
                    % store the 'sparse structure'
                    for si = 1:2
                        for sj = 1:2
                            obj.ii(index:index+8) = repmat((3*(obj.elem(i,si)-1)+1:3*obj.elem(i,si))',3,1);
                            obj.jj(index:index+8) = reshape(repmat(3*(obj.elem(i,sj)-1)+1:3*obj.elem(i,sj),3,1),9,1);
                            index = index + 9;
                        end
                    end
                end
            end
        end
        
        SetCurrentState(obj, x, v)
        
        SetMaterial(obj, k, Rho, elem, type)
        
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
        
        energy = totalEnergy(obj)
        
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
            for t = 1:obj.NS
                assert(obj.elemMaterialType(t) ~= 0); % check all elements have been initialized
            end
            obj.finalized = true;
        end
        function E = energy(obj)
            elastic_energy = 0;
            for i = 1:obj.NS
                Fs = obj.F(3*(i-1)+1:3*i,:);
                ks = obj.ks(i);
                elastic_energy = elastic_energy + 1/2 * ks * ((Fs'*Fs - 1)^2);
            end
            kinetic_energy = 0;
            for i = 1:obj.N
                vel = obj.v(3*(i-1)+1:3*i);
                kinetic_energy = kinetic_energy + 1/2 * 1 *(vel' * vel);
                % to do: change mass
            end
            potential = 0;
            for i = 1:obj.N
                h = obj.X(3*i) - obj.x(3*i);
                potential = potential - 9.8 * 1 * h;
                % to do: change mass
            end
            E = elastic_energy + potential + kinetic_energy;
        end
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

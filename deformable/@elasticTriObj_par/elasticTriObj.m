classdef elasticTriObj < elasticObj

    % elasticTriObj An elastic object represented by a triangle mesh
    %   Detailed explanation goes here
        
    properties (Constant)
        gravity = [0 -9.81]';
        
        G = [1 0; 0 1; -1 -1];
        Iv = eye(2);
        Im = eye(4);
        % commutation matrix
        Kmm = [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1];
        
        
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
        CGTri; % CG triangulation
        Dim = 2;
    end
    
    methods
        function obj = elasticTriObj(varargin)
            % Constructor of elasticTriObj
            %   Take in an undeformed tri mesh with
            %       nodeM: (undeformed) nodal position in material (#nodes
            %       by 2)
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
                    vol = det(T_node'*obj.G)/2; % undeformed volume from matrix determinant
                    assert(vol > 0, 'need to fix mesh orientation'); % if volume not positive, the input mesh need to be fixed
                    obj.W(i) = vol;
                    obj.T(4*(i-1)+1:4*i,:) = kron((obj.G * obj.DmINV(2*(i-1)+1:2*i,:))', obj.Iv); % definition: vec(F) = T * vec(x), or vec(dF) = T * vec(dx)
                    
                    for ti = 1:3
                        for tj = 1:3
                            obj.ii(index:index+3) = repmat((2*(obj.elem(i,ti)-1)+1:2*obj.elem(i,ti))',2,1);
                            obj.jj(index:index+3) = reshape(repmat(2*(obj.elem(i,tj)-1)+1:2*obj.elem(i,tj),2,1),4,1);
                            index = index + 4;
                        end
                    end
                    
                end
                
                obj.calculateGravity;
            end
            
        end
        
        f = ElasticForce(obj)
        
        function simple_vis(obj,ax,varargin)
            switch nargin
                case 1
                    triplot(obj.elem, obj.node(:,1),obj.node(:,2));
                case 2
                    axes(ax);
                    triplot(obj.elem, obj.node(:,1),obj.node(:,2));
                    obj.vis_handle = ax;
                case 3 % case where color is specified
                    color = varargin{1};
                    axes(ax);
                    hl = triplot(obj.elem, obj.node(:,1),obj.node(:,2),'Color',color);
                    obj.vis_handle = ax;     
            end
        end
        
        SetMaterial(obj, Y, P, Rho, type, a, b)
        
        ElasticForceDifferential(obj)
        
        K = StiffnessMatrix(obj)
        
        totalEnergy(obj)
        
        mesh_quality(obj)

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
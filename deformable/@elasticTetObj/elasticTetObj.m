classdef elasticTetObj < elasticObj
    %elasticTetObj An elastic object represented by a tet mesh
    %   Detailed explanation goes here
    
    properties (Constant)
        
        gravity = [0 0 -9.81]';
        
        
        G = [1 0 0; 0 1 0; 0 0 1; -1 -1 -1]; % element-wise operator that maps element nodal position in to edge matrix Ds or Dm
        Iv = eye(3);
        Im = eye(9);
        % commutation matrix
        Kmm =   ...
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
        Dim = 3;
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
                assert(size(input_nodeM,2) == obj.Dim);
                
                obj.N = size(input_nodeM,1);
                obj.NT = size(input_elem,1);
                
                obj.nodeM = input_nodeM;
                obj.elem = input_elem;
                
                obj.X = reshape(input_nodeM',obj.Dim*obj.N,1);
                obj.Dm = zeros(obj.Dim*obj.NT,obj.Dim);
                obj.DmINV = zeros(obj.Dim*obj.NT,obj.Dim);
                obj.F = zeros(obj.Dim*obj.NT,obj.Dim);
                obj.U = zeros(obj.Dim*obj.NT,obj.Dim);
                obj.V = zeros(obj.Dim*obj.NT,obj.Dim);
                obj.S = ones(obj.Dim*obj.NT,1);
                obj.R = zeros(obj.Dim*obj.NT,obj.Dim);
                obj.C = zeros(obj.Dim*obj.NT,obj.Dim);
                obj.normC = zeros(obj.NT,1);
                obj.M = sparse(size(obj.X,1),size(obj.X,1));
                obj.W = zeros(obj.NT,1);
                obj.T = zeros(9*obj.NT, 12);
                index = 1;
                load('local_k_ind.mat');
                for i = 1:obj.NT
                    T_node = input_nodeM(input_elem(i,:),:); % element nodal position in material space  (#nodes per elements by obj.Dim)
                    obj.Dm(obj.Dim*(i-1)+1:obj.Dim*i,:) = T_node'*obj.G;
                    obj.DmINV(obj.Dim*(i-1)+1:obj.Dim*i,:) = inv(T_node'*obj.G);
                    obj.W(i) =  abs(det(T_node'*obj.G)/6); % undeformed volume from matrix determinant
                    % using the negative determinant because of orientation
                    obj.T(9*(i-1)+1:9*i,:) = kron((obj.G * obj.DmINV(obj.Dim*(i-1)+1:obj.Dim*i,:))', obj.Iv); % definition: vec(F) = T * vec(x), or vec(dF) = T * vec(dx)
                    for e = input_elem(i,:)
                        for mi = (e-1)*obj.Dim+1:e*obj.Dim
                            obj.M(mi,mi) = obj.M(mi,mi)+obj.W(i)/4;
                        end
                    end
                    
                    for ti = 1:4
                        for tj = 1:4
                            obj.ii(index:index+8) = repmat((obj.Dim*(obj.elem(i,ti)-1)+1:obj.Dim*obj.elem(i,ti))',obj.Dim,1);
                            obj.jj(index:index+8) = reshape(repmat(obj.Dim*(obj.elem(i,tj)-1)+1:obj.Dim*obj.elem(i,tj),obj.Dim,1),9,1);
                            index = index + 9;
                        end
                    end
                    temp_ii = obj.ii(index-144:index-1);
                    temp_jj = obj.jj(index-144:index-1);
                    obj.ii(index-144:index-1) = temp_ii(local_k_ind);
                    obj.jj(index-144:index-1) = temp_jj(local_k_ind);
                end
                
                obj.calculateGravity;
            end
            
        end
        
        f = ElasticForce(obj)
        
        
        function simple_vis(obj,ax,varargin)
                    p = obj.node;
                    t = obj.elem;
                    bcol = [0.9,0.8,0.5];
            switch nargin
                case 1
                    tri1=surftri(p,t);
                    h=trimesh(tri1,p(:,1),p(:,2),p(:,3));
                    hold off
                    set(h,'facecolor',bcol,'edgecolor','k');
                case 2
                    axes(ax);

                    tri1=surftri(p,t);
                    h=trimesh(tri1,p(:,1),p(:,2),p(:,3));
                    hold off
                    set(h,'facecolor',bcol,'edgecolor','k');
                    obj.vis_handle = ax;
                case 3 % case where color is specified
                    bcol = varargin{1};
                    axes(ax);
                    
                    tri1=surftri(p,t);
                    h=trimesh(tri1,p(:,1),p(:,2),p(:,3));
                    hold off
                    set(h,'facecolor',bcol,'edgecolor','k');
                    obj.vis_handle = ax;     
            end
            obj.patch_handle = h;
        end
        
        [V_s, low_eig] = EigenFit_subspace(obj,K,M);
        [V_s, low_eig] = EigenFit_subspace2(obj,K,M);

        
        K = eigModification(obj,eig)
        
        [K,f] = eigModification_nonlinear_subspace(obj);

        
        f = eigModification_nonlinear_subspace_f_only(obj);

        [K,f] = eigModification_nonlinear_subspace_ratio_fix(obj);

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
        test_tet_orientation()
    end
end

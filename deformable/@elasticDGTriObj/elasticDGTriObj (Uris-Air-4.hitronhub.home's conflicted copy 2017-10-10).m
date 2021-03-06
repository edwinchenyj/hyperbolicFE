classdef elasticDGTriObj < handle

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
        CGTri; % CG triangulation
        
        DGN; % #DG nodes (3NT)
        DGnode; % DG nodal position in world (DGN by 2);
        DGnodeM; % DG nodal position in world (DGN by 2);
        DGelem; % DG elem list;
        DGM; % DG mass matrix
        DGK; % DG stiffnes matrix
        
        MapDGnode; % map from mesh nodal indices to indices in DGnode
        DGTri; % DG triangulation
        
        HalfEdge; % half edges
        Edge;
        MapHalfEdgeToElement;
        MapHalfEdge;
        DGHalfEdge; % half edges
        DGEdge;
        
        DGx;
        DGv;
        DGX;
        
        DGElement_ii;
        DGElement_jj;
        
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
        materialBlock = {}; % A cell array of indice to the elements of the same material type
        materialType = {}; % A cell array of type of material of each block
        materialBlockCount = 0; % number of different blocks of material
        elemMaterialType; % an array containing the material type of each element (NT by 1), 1 = neo-hookean, 2 = linear
        
        mu; % An array of mu to the elements 
        lambda; % An array of lambda to the elements 
        rho; % An array of rho to the elements 
        
        % TODO: need to get rid of this
        Y;
        P;
        Rho;
        
        % for ghost sub-points
        sub_objs;
        
        % handle to the axis for visualization
        vis_handle;
        
        % indices for fixed points (for ghost points)
        ind_fix;
        ind_apex;
        ind_apex2;
        ind_apex3;
    end
    
    methods
        function obj = elasticDGTriObj(varargin)
            % Constructor 
            %   Take in an undeformed tri mesh with
            %       nodeM: (undeformed) nodal position in material (#nodes
            %       by 3)
            %       elem: element label (#element by 3) [N1 N2 N3]
            if nargin == 1 && isa(varargin{1},'elasticDGTriObj')
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
                
                obj.DGelem = zeros(size(obj.elem));
                obj.DGN = 3 * obj.NT; % #DG nodes (3NT)
                obj.DGnode = zeros(obj.DGN,2); % DG nodal position in world (DGN by 2);
                obj.DGnodeM = zeros(obj.DGN,2); % DG nodal position in world (DGN by 2);
                
                obj.MapDGnode = cell(obj.N,1);
                
                % initialize the vector to zeros by flatten DGnode
                obj.DGx = obj.DGnode(:);
                obj.DGv = obj.DGnode(:);
                obj.DGX = obj.DGnode(:);
                
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
                obj.DGM = sparse(size(obj.DGX,1),size(obj.DGX,1));
                obj.W = zeros(obj.NT,1);
                obj.T = zeros(4*obj.NT, 6);
                
                obj.Edge = []; % don't know the edge size when the mesh is created
                obj.DGEdge = [];
                obj.MapHalfEdgeToElement = [];
                obj.MapHalfEdge = [];
                
                
                obj.CGTri = triangulation(input_elem,input_nodeM);
                
                index = 1;
                edge_index = 1;
                
                for i = 1:obj.NT
                    T_node = input_nodeM(input_elem(i,:),:); % element nodal position in material space  (#nodes per elements by 2)
                    for iNode = 1:3
                        if iNode == 1
                            obj.MapDGnode{input_elem(i,iNode)} = [obj.MapDGnode{input_elem(i,iNode)} 3*(i-1)+1]; 
                        elseif iNode == 2
                            obj.MapDGnode{input_elem(i,iNode)} = [obj.MapDGnode{input_elem(i,iNode)} 3*(i-1)+2];
                        else
                            obj.MapDGnode{input_elem(i,iNode)} = [obj.MapDGnode{input_elem(i,iNode)} 3*i]; 
                        end
                    end
                    obj.DGnodeM(3*(i-1)+1:3*i,:) = T_node;
                    obj.DGelem(i,:) = (3*(i-1)+1:3*i);
                    obj.HalfEdge = [obj.Edge; intput_elem(i,1) intput_elem(i,2); intput_elem(i,2) intput_elem(i,3); intput_elem(i,3) intput_elem(i,1)];
                    obj.DGEdge = [obj.DGEdge; 3*(i-1)+1 3*(i-1)+2; 3*(i-1)+2 3*i; 3*i 3*(i-1)+1];
                    
                    % the map from half edges to the element it belongs to
                    obj.MapHalfEdgeToElement = [obj.MapHalfEdgeToElement;...
                        edge_index i;...
                        edge_index + 1 i;...
                        edge_index + 2 i];
                    
                    k1 = find(sum(obj.HalfEdge == [intput_elem(i,2) intput_elem(i,1)],2) == 2);
                    if (~isempty(k1))
                        obj.HalfEdge(k1,:) = [obj.HalfEdge(k1,:) intput_elem(i,2) intput_elem(i,1)];
                        obj.HalfEdge(edge_index,:) = [obj.HalfEdge(k1,:) intput_elem(i,2) intput_elem(i,1)];

                        obj.DGEdge(k1,:) = [obj.DGEdge(k1,:) 3*(i-1)+1 3*(i-1)+2];
                    end
                    k2 = find(sum(obj.HalfEdge == [intput_elem(i,3) intput_elem(i,2)],2) == 2);
                    if (~isempty(k2))
                        obj.HalfEdge(k2,:) = [obj.HalfEdge(k2,:) intput_elem(i,3) intput_elem(i,2)];
                    end
                    k3 = find(sum(obj.HalfEdge == [intput_elem(i,1) intput_elem(i,3)],2) == 2);
                    if (~isempty(k3))
                        obj.HalfEdge(k3,:) = [obj.HalfEdge(k3,:) intput_elem(i,1) intput_elem(i,3)];
                    end
                    
                    edge_index = edge_index + 3;
                    
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
                            obj.DGElement_ii(index:index+3) = repmat((2*(obj.DGelem(i,ti)-1)+1:2*obj.DGelem(i,ti))',2,1);
                            obj.DGElement_jj(index:index+3) = reshape(repmat(2*(obj.DGelem(i,tj)-1)+1:2*obj.DGelem(i,tj),2,1),4,1);
                            index = index + 4;
                        end
                    end
                end
                
                
                obj.DGTri = triangulation(obj.DGelem, obj.DGnodeM);
                obj.DGX = obj.CGxToDGx(obj.X);
            end
            
        end
        
        SetCurrentState(obj, x)
        SetCurrentDGState(obj, DGx)
        
        SetMaterial(obj, Y, P, Rho, elem, type)
        
        function DGx = CGxToDGx(obj,CGx)
            DGx = zeros(size(obj.DGx));
            for i = 1:obj.N
                Nx = CGx(obj.Dim*(i-1)+1:obj.Dim*i);
                for iDG = obj.MapDGnode{i}
                    DGx(obj.Dim*(iDG-1)+1: obj.Dim*iDG) = Nx;
                end
            end
        end
        
%         function DGDx = CGDxToDGDx(obj,CGDx)
%         end
        
        function ha = init_vis(obj)
            ha = axes;
            obj.vis_handle = ha;
            hold(obj.vis_handle,'on');
        end
        
        function simple_vis(obj,ax)
            axes(ax);
            triplot(obj.elem,obj.nodeM(:,1),obj.nodeM(:,2));
            obj.vis_handle = ax;
%             ha = ax;
        end
        function current_vis(obj,ax)
            axes(ax);
            triplot(obj.elem,obj.node(:,1),obj.node(:,2));
            obj.vis_handle = ax;
%             ha = ax;
        end
        
        function simple_vis_sub(obj,ax)
            for t = 1:obj.NT
                local_elem = obj.sub_objs(t).elem;
                local_nodeM = obj.sub_objs(t).nodeM;
                
                axes(ax);
                triplot(local_elem,local_nodeM(:,1),local_nodeM(:,2));
            end
            
        end
        
        function DG_vis(obj,ax)
            axes(ax);
            triplot(obj.DGelem, obj.DGnodeM(:,1),obj.DGnodeM(:,2));
            obj.vis_handle = ax;
        end
        
                
        function DG_current_vis(obj,ax)
            axes(ax);
            triplot(obj.DGelem, obj.DGnode(:,1),obj.DGnode(:,2));
            obj.vis_handle = ax;
        end
        
        
        
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
        
        f = subElasticForce(obj)
        
        f = ElasticForce(obj)
        
        function out = ElasticForceWFixDx(obj,dx_free)
            Dx = obj.x - obj.X;
            Dx(~obj.ind_fix) = dx_free;
            obj.SetCurrentState(Dx);
            force = obj.ElasticForce;
            out = force(~obj.ind_fix);
        end
        
        df = ElasticForceDifferential(obj, dx)
        
        K = StiffnessMatrix(obj)
        
        energy = totalEnergy(obj)
        
        initGhostPoints(obj)
        
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
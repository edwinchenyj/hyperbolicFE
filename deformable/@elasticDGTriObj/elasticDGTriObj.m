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
        K0 =[];
        CGTri; % CG triangulation
        
        DGN; % #DG nodes (3NT)
        DGnode; % DG nodal position in world (DGN by 2);
        DGnodeM; % DG nodal position in world (DGN by 2);
        DGelem; % DG elem list;
        DGM; % DG mass matrix
        DGK; % DG stiffnes matrix
        DGK0 = [];
        MapDGnode; % map from mesh nodal indices to indices in DGnode
        DGTri; % DG triangulation
        
        HalfEdge; % half edges
        Edge;
        MapHalfEdgeToElement;
        MapHalfEdge;
        DGHalfEdge; % half edges
        HalfEdgeLength;
        DGEdge;
        DGInterface; % DG interface
        
        DGx;
        DGv;
        DGX;
        DGb;
        
        eta = 1e-1;
        
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
                obj.HalfEdgeLength = zeros(0,1);
                obj.HalfEdge = zeros(0,4); % don't know the edge size when the mesh is created
                obj.DGEdge = zeros(0,4);
                
                obj.MapHalfEdgeToElement = zeros(0,1);
                obj.MapHalfEdge = zeros(0,1);
                
                % map from cg nodes to dg nodes
                obj.MapDGnode = cell(obj.N,1);
                
                obj.DGInterface = [];
                
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
                    obj.HalfEdge = [obj.HalfEdge; input_elem(i,1) input_elem(i,2) 0 0; input_elem(i,2) input_elem(i,3) 0 0; input_elem(i,3) input_elem(i,1) 0 0];
                    obj.DGEdge = [obj.DGEdge; 3*(i-1)+1 3*(i-1)+2 0 0; 3*(i-1)+2 3*i 0 0; 3*i 3*(i-1)+1 0 0 ];
                    obj.HalfEdgeLength = [obj.HalfEdgeLength; norm(input_nodeM(input_elem(i,1),:)-input_nodeM(input_elem(i,2),:));...
                        norm(input_nodeM(input_elem(i,2),:)-input_nodeM(input_elem(i,3),:));...
                        norm(input_nodeM(input_elem(i,3),:)-input_nodeM(input_elem(i,1),:))];
                    % the map from half edges to the element it belongs to
                    obj.MapHalfEdgeToElement = [obj.MapHalfEdgeToElement;
                        i;...
                        i;...
                        i];
                    
                    % find if the same edge has been added as the
                    % neighboring half edge
                    %                     k1 = find(sum(obj.HalfEdge(:,1:2) == [input_elem(i,2) input_elem(i,1)],2) == 2);
                    k1 = find(ismember(obj.HalfEdge(:,1:2), [input_elem(i,2) input_elem(i,1)],'rows'));
                    if (~isempty(k1))
                        obj.HalfEdge(k1,3:4) = [input_elem(i,2) input_elem(i,1)];
                        obj.HalfEdge(edge_index,3:4) = obj.HalfEdge(k1,1:2);
                        obj.DGEdge(k1,3:4) = [input_elem(i,2) input_elem(i,1)];
                        obj.DGEdge(edge_index,3:4) = obj.DGEdge(k1,1:2);
                        obj.MapHalfEdge(k1) = edge_index;
                        obj.MapHalfEdge(edge_index) = k1;
                    end
                    k2 = find(ismember(obj.HalfEdge(:,1:2), [input_elem(i,3) input_elem(i,2)],'rows'));
                    if (~isempty(k2))
                        obj.HalfEdge(k2,3:4) = [input_elem(i,3) input_elem(i,2)];
                        obj.HalfEdge(edge_index+1,3:4) = obj.HalfEdge(k2,1:2);
                        obj.DGEdge(k2,3:4) = [input_elem(i,3) input_elem(i,2)];
                        obj.DGEdge(edge_index+1,3:4) = obj.DGEdge(k2,1:2);
                        obj.MapHalfEdge(k2) = edge_index + 1;
                        obj.MapHalfEdge(edge_index + 1) = k2;
                    end
                    k3 = find(ismember(obj.HalfEdge(:,1:2), [input_elem(i,1) input_elem(i,3)],'rows'));
                    if (~isempty(k3))
                        obj.HalfEdge(k3,3:4) = [input_elem(i,1) input_elem(i,3)];
                        obj.HalfEdge(edge_index+2,3:4) = obj.HalfEdge(k3,1:2);
                        obj.DGEdge(k3,3:4) = [input_elem(i,1) input_elem(i,3)];
                        obj.DGEdge(edge_index+2,3:4) = obj.DGEdge(k3,1:2);
                        obj.MapHalfEdge(k3) = edge_index + 2;
                        obj.MapHalfEdge(edge_index + 2) = k3;
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
                
                %                 for e = 1:size(obj.HalfEdge,1)
                %                     if obj.HalfEdge(e,3) ~= 0
                %                         obj.DGInterface =  [obj.DGInterface];
                %                     end
                %                 end
                
                
                
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
        
        
        function DG_current_vis(obj,ax,varargin)
            switch nargin
                case 1
                    triplot(obj.DGelem, obj.DGnode(:,1),obj.DGnode(:,2));
                case 2
                    axes(ax);
                    triplot(obj.DGelem, obj.DGnode(:,1),obj.DGnode(:,2));
                    obj.vis_handle = ax;
                case 3 % case where color is specified
                    color = varargin{1};
                    axes(ax);
                    hl = triplot(obj.DGelem, obj.DGnode(:,1),obj.DGnode(:,2),'Color',color);
                    obj.vis_handle = ax;
                    
            end
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
        
        f = DGElasticForce(obj)
        f = DGInterfaceElasticForce(obj);
        
        DGCalculateConst(obj)
        
        function out = ElasticForceWFixDx(obj,dx_free)
            Dx = obj.x - obj.X;
            Dx(~obj.ind_fix) = dx_free;
            obj.SetCurrentState(Dx);
            force = obj.ElasticForce;
            out = force(~obj.ind_fix);
        end
        
        df = ElasticForceDifferential(obj, dx)
        
        K = StiffnessMatrix(obj)
        
        K = DGElmentStiffnessMatrix(obj)
        
        K = DGStiffnessMatrix(obj)
        
        K = DGInterfaceStiffnessMatrix(obj)
        
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
        
        function out = quadratic_line_int(p0,p1,A)
            len = norm(p1 - p0);
            x0 = p0(1); y0 = p0(2);
            x1 = p1(1); y1 = p1(2);
            a11 = A(1,1); a12 = A(1,2); a21 = A(2,1); a22 = A(2,2);
            out = len * ( a11/3 * (x1*x1 + x1*x0 + x0*x0)...
                + a22/3 * (y1*y1 + y1*y0 + y0*y0)...
                + (a12 + a21) * (1/6*x0*y1 + 1/6*y0*x1 + 1/3*x1*y1 + 1/3*x0*y0));
        end
        
        function out = linear_line_int(p0,p1,b)
            len = norm(p1 - p0);
            x0 = p0(1); y0 = p0(2);
            x1 = p1(1); y1 = p1(2);
            b1 = b(1); b2 = b(2);
            out = len * (b1*x0 + b2*y0 + 1/2*b1*(x1-x0) + 1/2*b2*(y1-y0));
        end
        
        function out = const_line_int(p0,p1,c)
            len = norm(p1 - p0);
            out = c*len;
        end
    end
    
end
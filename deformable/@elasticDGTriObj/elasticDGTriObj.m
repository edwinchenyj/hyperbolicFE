classdef elasticDGTriObj < elasticObj
    
    % elasticTriObj An elastic object represented by a triangle mesh
    %   Detailed explanation goes here
    
    properties (Constant)
        
        gravity = [0 -9.81];
        
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
        Dim = 2;
        
        CGTri; % CG triangulation
        CGN;
        
        DGN; % #DG nodes (3NT)
        DGnode; % DG nodal position in world (DGN by 2);
        DGnodeM; % DG nodal position in world (DGN by 2);
        DGelem; % DG elem list;
        DGM; % DG mass matrix
        DGK; % DG stiffnes matrix
        DGK0 = [];
        DGKi = [];
        MapDGnode; % map from mesh nodal indices to indices in DGnode
        DGTri; % DG triangulation
        
        DGBZ;
        DGIP;
        
        HalfEdge; % half edges
        Edge;
        MapHalfEdgeToElement;
        MapHalfEdge;
        DGHalfEdge; % half edges
        HalfEdgeLength;
        HalfEdgeUnitOutNormal;
        DGEdge;
        DGInterface; % DG interface
        
        DGx;
        DGv;
        DGX;
        DGb;
        
        eta = 1e-1;
        
        DGElement_ii;
        DGElement_jj;
        
        InterfacePhi = [];
        InterfaceExp = [];
        
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
                obj.CGN = size(input_nodeM,1); % number of CG nodes
                

                obj.NT = size(input_elem,1);
                obj.N = 3 * obj.NT; % number of DG nodes
                obj.elem = input_elem;
                
%                 obj.DGelem = zeros(size(input_elem));
                obj.node = zeros(obj.N,2); % DG nodal position in world (DGN by 2);
                obj.nodeM = zeros(obj.N,2); % DG nodal position in world (DGN by 2);
                
                % initialize the vector to zeros by flatten DGnode
                obj.x = obj.node(:);
                obj.v = obj.node(:);
                obj.X = obj.node(:);
                
                % initialize DG to be BZ and IP
                obj.DGBZ = true;
                obj.DGIP = true;
                
                CGX = reshape(input_nodeM',2*obj.CGN,1); % rest state X
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
%                 obj.DGM = sparse(size(obj.DGX,1),size(obj.DGX,1));
                obj.W = zeros(obj.NT,1);
                obj.T = zeros(4*obj.NT, 6);
                obj.HalfEdgeLength = zeros(0,1);
                obj.HalfEdgeUnitOutNormal = zeros(0,1);
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
                    vol = det(T_node'*obj.G)/2; % undeformed volume from matrix determinant
                    assert(vol > 0, 'need to fix mesh orientation'); % if volume not positive, the input mesh need to be fixed
                    obj.W(i) =  vol;
                    obj.nodeM(3*(i-1)+1:3*i,:) = T_node;
                    obj.elem(i,:) = (3*(i-1)+1:3*i);
                    obj.HalfEdge = [obj.HalfEdge; input_elem(i,1) input_elem(i,2) 0 0; input_elem(i,2) input_elem(i,3) 0 0; input_elem(i,3) input_elem(i,1) 0 0];
                    obj.DGEdge = [obj.DGEdge; 3*(i-1)+1 3*(i-1)+2 0 0; 3*(i-1)+2 3*i 0 0; 3*i 3*(i-1)+1 0 0 ];
                    dE1 = input_nodeM(input_elem(i,2),:)-input_nodeM(input_elem(i,1),:);
                    dE2 = input_nodeM(input_elem(i,3),:)-input_nodeM(input_elem(i,2),:);
                    dE3 = input_nodeM(input_elem(i,1),:)-input_nodeM(input_elem(i,3),:);
                    obj.HalfEdgeLength = [obj.HalfEdgeLength; norm(dE1);...
                        norm(dE2);...
                        norm(dE3)];
                    obj.HalfEdgeUnitOutNormal = [obj.HalfEdgeUnitOutNormal; [dE1(2) -dE1(1)]/ norm([dE1(2) -dE1(1)]);...
                        [dE2(2) -dE2(1)]/norm([dE2(2) -dE2(1)]);...
                        [dE3(2) -dE3(1)]/norm([dE3(2) -dE3(1)])];
                    % the map from half edges to the element it belongs to
                    obj.MapHalfEdgeToElement = [obj.MapHalfEdgeToElement;
                        i;...
                        i;...
                        i];
                    
                    % find if the same edge has been added as the
                    % neighboring half edge 
                    % search in the reverse order because of orientation
                    %                     k1 = find(sum(obj.HalfEdge(:,1:2) == [input_elem(i,2) input_elem(i,1)],2) == 2);
                    k1 = find(ismember(obj.HalfEdge(:,1:2), [input_elem(i,2) input_elem(i,1)],'rows'));
                    if (~isempty(k1))
                        obj.HalfEdge(k1,3:4) = [input_elem(i,1) input_elem(i,2)];
                        obj.HalfEdge(edge_index,3:4) = obj.HalfEdge(k1,1:2);
                        obj.DGEdge(k1,3:4) = [input_elem(i,1) input_elem(i,2)];
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

                    obj.T(4*(i-1)+1:4*i,:) = kron((obj.G * obj.DmINV(2*(i-1)+1:2*i,:))', obj.Iv); % definition: vec(F) = T * vec(x), or vec(dF) = T * vec(dx)
                    
                    for ti = 1:3
                        for tj = 1:3
                            obj.ii(index:index+3) = repmat((2*(input_elem(i,ti)-1)+1:2*input_elem(i,ti))',2,1);
                            obj.jj(index:index+3) = reshape(repmat(2*(input_elem(i,tj)-1)+1:2*input_elem(i,tj),2,1),4,1);
                            obj.DGElement_ii(index:index+3) = repmat((2*(obj.elem(i,ti)-1)+1:2*obj.elem(i,ti))',2,1);
                            obj.DGElement_jj(index:index+3) = reshape(repmat(2*(obj.elem(i,tj)-1)+1:2*obj.elem(i,tj),2,1),4,1);
                            index = index + 4;
                        end
                    end
                end
                
                obj.DGTri = triangulation(obj.elem, obj.nodeM);
                obj.X = obj.CGxToDGx(CGX);
                obj.calculateGravity;
            end
            
        end
        
        SetCurrentState(obj, x)
        
        SetMaterial(obj, Y, P, Rho, type, a, b)
        
        function DGx = CGxToDGx(obj,CGx)
            DGx = zeros(size(obj.x));
            for i = 1:obj.CGN
                Nx = CGx(obj.Dim*(i-1)+1:obj.Dim*i);
                for iDG = obj.MapDGnode{i}
                    DGx(obj.Dim*(iDG-1)+1: obj.Dim*iDG) = Nx;
                end
            end
        end
        
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
        
        phi = DGEXPINTInterfacePhi(obj,dt);
        
        f = ElasticForce(obj)
        
        f = DGInterfaceElasticForce(obj);
        
        DGCalculateConst(obj)
        
        df = ElasticForceDifferential(obj, dx)
        
        K = StiffnessMatrix(obj)
        
        K = DGElmentStiffnessMatrix(obj)
        
        K = DGInterfaceStiffnessMatrix(obj)
        
        energy = totalEnergy(obj)
        
        
        mesh_quality(obj)

        P = Stress(obj, t)
        C = FourtOrderTensor(obj, t)
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
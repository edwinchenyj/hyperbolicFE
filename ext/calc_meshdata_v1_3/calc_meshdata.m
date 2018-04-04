function md = calc_meshdata(dim,p,ef,t)

%LOAD_MESHDATA   Calculates stuff needed for integration in 2D/3D.
%
% SYNTAX (2D):   md = load_meshdata(p,e,t)
%
% IN:   dim         the dimension of the mesh
%       p,e,t       the variables given by the initmesh/poimesh
%                   functions found in the PDEtool toolbox
%
% SYNTAX (3D):   md = load_meshdata(p,f,t)
%
% IN:   dim         the dimension of the mesh
%       p,f,t       the variables given by the read_tetgenmesh
%                   function
%
% OUT:  md.         the meshdata structure which contains:
%        (first the input variables)
%        p          points                   
%        e          boundary edges by points (only in 2D)
%        f          boundary edges by points (only in 3D)
%        t          elements by points
%
%        (then the calculated variables)
%        B_K        the affine mapping from the unit triangle:
%        b_K            F_K(x) = B_K * x + b_K
%        B_K_inv    the inverse of B_K
%        B_K_det    the determinant of B_K
%
%        areas      the areas of triangles (only in 2D)
%        vols       the areas of triangles (only in 3D)
%

%
%   Author:  Immanuel Anjam, University of Jyvaskyla
%   Contact: immanuel.anjam@gmail.com
%   Date:    26.09.2014
%   Version: 1.2
%

if ( dim == 2 )

    md.p = p;                   % gather the input
    md.e = ef;
    md.t = t;

    nelems = size(md.t,2);      % initialize
    md.B_K = zeros(2,2,nelems);

    A = md.p(1:2,md.t(1,:));    % points defining the triangle
    B = md.p(1:2,md.t(2,:));
    C = md.p(1:2,md.t(3,:));
    a = B - A;                  % vectors defining the triangle
    b = C - A;

    md.B_K(:,1,:) = a;          % the affine mapping
    md.B_K(:,2,:) = b;
    md.b_K        = A;

    md.B_K_det = a(1,:).*b(2,:) - a(2,:).*b(1,:);       % determinant
    const      = md.B_K_det.^(-1);
    md.B_K_inv(1,1,:) =   b(2,:) .* const;              % inverse of B_K
    md.B_K_inv(1,2,:) = - b(1,:) .* const;
    md.B_K_inv(2,1,:) = - a(2,:) .* const;
    md.B_K_inv(2,2,:) =   a(1,:) .* const;

    md.areas          = 1/2 * abs( md.B_K_det );        % areas
    
elseif ( dim == 3 )
    
    md.p = p;                   % gather the input
    md.f = ef;
    md.t = t;

    nelems = size(md.t,2);      % initialize
    md.B_K = zeros(3,3,nelems);
    
    A = md.p(1:3,md.t(1,:));    % points defining the triangle
    B = md.p(1:3,md.t(2,:));
    C = md.p(1:3,md.t(3,:));
    D = md.p(1:3,md.t(4,:));
    a = B - A;                  % vectors defining the triangle
    b = C - A;
    c = D - A;

    md.B_K(:,1,:) = a;          % the affine mapping
    md.B_K(:,2,:) = b;
    md.B_K(:,3,:) = c;
    md.b_K = A;

    cp = [ b(2,:).*c(3,:) - b(3,:).*c(2,:) ;    % crossproduct b x c
           b(3,:).*c(1,:) - b(1,:).*c(3,:) ;
           b(1,:).*c(2,:) - b(2,:).*c(1,:) ];
    md.B_K_det = sum( a .* cp );                % determinant

    const = md.B_K_det.^(-1);                   % calculate the inverse matrix
    md.B_K_inv(1,1,:) = ( b(2,:).*c(3,:) - b(3,:).*c(2,:) ) .* const;
    md.B_K_inv(1,2,:) = ( c(1,:).*b(3,:) - c(3,:).*b(1,:) ) .* const;
    md.B_K_inv(1,3,:) = ( b(1,:).*c(2,:) - b(2,:).*c(1,:) ) .* const;
    md.B_K_inv(2,1,:) = ( c(2,:).*a(3,:) - c(3,:).*a(2,:) ) .* const;
    md.B_K_inv(2,2,:) = ( a(1,:).*c(3,:) - a(3,:).*c(1,:) ) .* const;
    md.B_K_inv(2,3,:) = ( c(1,:).*a(2,:) - c(2,:).*a(1,:) ) .* const;
    md.B_K_inv(3,1,:) = ( a(2,:).*b(3,:) - a(3,:).*b(2,:) ) .* const;
    md.B_K_inv(3,2,:) = ( b(1,:).*a(3,:) - b(3,:).*a(1,:) ) .* const;
    md.B_K_inv(3,3,:) = ( a(1,:).*b(2,:) - a(2,:).*b(1,:) ) .* const;

    md.vols = 1/6 * abs( md.B_K_det );                 % volumes
    
else
    
    error('calc_meshdata: Only dimensions 2 and 3, please.')
    
end
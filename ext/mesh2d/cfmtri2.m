function [vert,econ,tria] = cfmtri2(vert,econ)
%CFMTRI2 compute a conforming 2-simplex Delaunay triangulat-
%ion in the two-dimensional plane.
%   [VERT,CONN,TRIA]=CFMTRI2(VERT,CONN) computes the confor-
%   ming Delaunay trianguation, given the points VERT, and
%   edge constraints CONN. New points are inserted to bisect
%   edge constraints until all are recovered. VERT is a 
%   V-by-2 array of XY coordinates to be triangulated, TRIA 
%   is a T-by-3 array of vertex indexing, with each row 
%   defining a triangle, such that VERT(TRIA(II,1),:), 
%   VERT(TRIA(II,2),:) and VERT(TRIA(II,3),:) are the coord-
%   inates of the II-TH triangle. CONN is a C-by-2 array of
%   constraining edges, where each row defines an edge, as
%   per the TRIA array.
%
%   See also DELTRI2, DELAUNAYN

%   Darren Engwirda : 2017 --
%   Email           : engwirda@mit.edu
%   Last updated    : 24/03/2017

%---------------------------------------------- basic checks    
    if ( ~isnumeric(vert) || ...
         ~isnumeric(econ) )
        error('cfmtri2:incorrectInputClass' , ...
            'Incorrect input class.') ;
    end

%---------------------------------------------- basic checks
    if (ndims(vert) ~= +2 || ndims(econ) ~= +2)
        error('cfmtri2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end
    if (size(vert,2)~= +2 || size(econ,2)~= +2)
        error('cfmtri2:incorrectDimensions' , ...
            'Incorrect input dimensions.');
    end
    
%-- the DELAUNAYN routine is *not* well-behaved numerically,
%-- so explicitly re-scale the problem about [-1,-1; +1,+1]. 
    vmax = max(vert,[],1) ;
    vmin = min(vert,[],1) ;
    
    vdel = vmax - vmin;
    vdel = mean(vdel) ;
    vdel = vdel * +.5 ;
    
    vmid = vmax + vmin;
    vmid = vmid * +.5 ;
    
    vert = vert - vmid;
    vert = vert / vdel;
    
%-- keep bisecting edge constraints until they are all reco-
%-- vered! 
    while (true)

    %----------------- un-constrained delaunay triangulation
        tria = delaunay( ...
            vert(:,1),vert(:,2)) ;

        nv = size(vert,+1);
        nt = size(tria,+1);     
        
    %----------------------------- build non-unique edge-set
        ee = zeros(nt*3,2);
        ee((1:nt)+nt*0,:) = tria(:,[1,2]);
        ee((1:nt)+nt*1,:) = tria(:,[2,3]);
        ee((1:nt)+nt*2,:) = tria(:,[3,1]);

    %----------------- find constraints within tria-edge set
       [in] = setset2(econ,ee) ;

    %----------------------------- done when have contraints
        if (all(in)), break; end

    %----------------------------- un-recovered edge centres
        vm = vert(econ(~in,1),:) ...
           + vert(econ(~in,2),:) ;
        vm = vm * +.5 ;
        
    %----------------------------- un-recovered edge indexes
        ev = nv+(1:size(vm,1))';
        en = [econ(~in,+1), ev;
              econ(~in,+2), ev];
        
    %----------------------------- push new vert/edge arrays    
        vert = [vert( :,:); vm];
        econ = [econ(in,:); en];
 
    end

%--------------------------------- undo geomertic re-scaling
    vert = vert * vdel ;
    vert = vert + vmid ;

end




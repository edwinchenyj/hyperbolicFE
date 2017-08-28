%% control the qualities of 2d tri mesh of a right angle triangle
regenerate = false;


% type_list = {'triangle' 'circle' 'rectangleCircularHole' 'rect','poly','square' 'ellipse'};
type_list = {'triangle'};

% edge length list
el_list = [0.2 0.1 0.05 0.02];
% max area list
maxA_list = [1  0.3 0.2 0.1 0.05 0.01]';
% maxA_list = [2]
% now each area for quality control in triangle corresponds to each edge length from distmesh
nodeM = [0 0; 1 0; 0 1];
elem = [1 2 3];

type = 'rightangletriangle';
for i = 1:length(maxA_list)
    maxA = maxA_list(i);
    
    meshname = sprintf('simData%ctrimeshes%c%s maxA%.e tritool%s',fs,fs,type, maxA);
    if or(regenerate,exist([meshname '.mat'], 'file') ~= 2)
        
        [TV,TF,TN] = triangle(nodeM, elem, 'MaxArea',maxA,'Quality',30);
        
        p = TV;
        t = TF;
        [p,t] = fixmesh(p,t);
        % fstats(p,t);
        
        nodeM = p;
        elem = t;
        
    else
        load([meshname '.mat'], 'nodeM', 'elem');
    end
    colormap jet
    % figure
    trimesh(elem,nodeM(:,1),nodeM(:,2));
    axis equal
    grid on
    % figure
    % simpplot(nodeM,elem)
    %         pause
    save([meshname '.mat'],'nodeM','elem')
    
    print(meshname,'-dpng')
end


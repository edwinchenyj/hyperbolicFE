%% control the qualities of 2d tri mesh from distmesh using triangle with gptoolbox
regenerate = false;

fs = filesep;

type_list = {'triangle' 'circle' 'rectangleCircularHole' 'rect','poly','square' 'ellipse'};
type_list = {'triangle'};

% edge length list
el_list = [0.2 0.1 0.05 0.02];
% max area list
% maxA_list = [0.2 0.05 0.03 0.02]'/50;
maxA_list = [0.1]
% now each area for quality control in triangle corresponds to each edge length from distmesh


for i_type = 1:length(type_list)
    type = type_list{i_type};
    for i = 1:length(maxA_list)
        el = el_list(i);
        maxA = maxA_list(i);
        fprintf('%s\n',type);
        oldmeshname = sprintf('simData%ctrimeshes%c%s el%.e',fs,fs,type, el);
        
        meshname = sprintf('simData%ctrimeshes%c%s maxA%.e%s',fs,fs,type, maxA,'triangle');
        if or(regenerate,exist([meshname '.mat'], 'file') ~= 2)
            
            load([oldmeshname '.mat'], 'nodeM', 'elem');
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
        s_alpha = sprintf('%f',el);
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
    
    
end

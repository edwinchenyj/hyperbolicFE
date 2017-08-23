clear all


close all


close all
clear all

draw = true;
rerun = true;

type_list = {'triangle' 'circle' 'rectangleCircularHole' 'rect','poly','square' 'ellipse'};
% type_list = {'ellipse'}
% type_list = {'rectangleCircularHole'};
type_list = {'poly'};

% max area list
maxA_list = [0.2 0.05 0.03 0.02]'/50;

maxA_list = [0.03]'/50;

Xsource_indices = [1];
Ysource_indices = [1];

fs = filesep;
mass_flag = true;
for i_type = 1:length(type_list)
    type = type_list{i_type};
    
    
    total_mass_list = zeros(length(maxA_list),1);
    
    for i_maxA = 1:length(maxA_list)
        maxA = maxA_list(i_maxA);
        fprintf('%s\n',type);
        meshname = sprintf('simData%ctrimeshes%c%s maxA%.e tritool%s',fs,fs,type, maxA,'triangle');
        
        if exist([meshname '.mat'], 'file') ~= 2
            disp('mesh does not exist')
            break
        else
            load([meshname '.mat'], 'nodeM', 'elem');
            
        end
        %
        % % check tet orientation
        % for i = 1:size(elem,1)
        %     T = [nodeM(elem(i,:),:), ones(4,1)];
        %     assert(det(T)<0)
        % end
        
%         temp = nodeM(:,3);
%         nodeM(:,3) = nodeM(:,2);
%         nodeM(:,2) = temp;
%         nodeM = nodeM';
%         elem = elem';
        
        % tree = opcodemesh(bunny_vert,bunny_face);
        trimesh(elem,nodeM(:,1),nodeM(:,2));
        hold on
        hp = plot(nodeM(1,1),nodeM(2,1),'r.','MarkerSize',10);
        
        % n = 5;
        ht = text(nodeM(1,1),nodeM(2,1),num2str(1),'FontSize',8);
        %
        % % for i = 1:size(v,2)
        
        
        for j = 0:6
            n = 100*j+1;
            for i = n:n+100
                % for i = [256 2465]
                % for i = [256 964 966 965]
                %     a.String = num2str(i);
                %     a.Position = v(:,i);
                
                disp(i)
                disp(nodeM(i,:))
                text(nodeM(i,1),nodeM(i,2),num2str(i),'FontSize',8);
                %     hp.XData = v(1,i);
                %     hp.YData = v(2,i);
                %     hp.ZData = v(3,i);
                plot(nodeM(i,1),nodeM(i,2),'r.','MarkerSize',10)
                
                %     drawnow
                %     pause(0.1)
            end
            pause
        end
        %
    end
end
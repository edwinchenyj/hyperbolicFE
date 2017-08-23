%% generate 2d tri mesh from distmesh
regenerate = false;


type_list = {'triangle' 'circle' 'rectangleCircularHole' 'rect','poly','square' 'ellipse'};
type_list = {'triangle'};

% edge length list
el_list = [0.1 0.05 0.025 0.02 0.015];
el_list = [0.2 0.1 0.05 0.02];
el_list = [0.005]

% everything in a fixed bounding box
boundingBox = [-0 -0; 1 1];

for i_type = 1:length(type_list)
    type = type_list{i_type};
    for i = 1:length(el_list)
        el = el_list(i);
        fprintf('%s\n',type);
        meshname = sprintf('simData%ctrimeshes%c%s el%.e',fs,fs,type, el);
        
        if or(regenerate,exist([meshname '.mat'], 'file') ~= 2)
            rng('default'); % Always the same results
            rng(1);
            
            if strcmp(type,'triangle')
                pv=[0 0; 1 0; 0 1; 0 0];
                [p,t]=distmesh2d(@dpoly,@huniform,el,boundingBox,pv,pv);
            elseif strcmp(type,'square')
                pv=[0 0; 1 0; 1 1; 0 1; 0 0];
                [p,t]=distmesh2d(@dpoly,@huniform,el,boundingBox,pv,pv);
            elseif strcmp(type,'rect')
                pv=[0 0; 3 0; 3 1; 0 1; 0 0]/3;
                [p,t]=distmesh2d(@dpoly,@huniform,el,boundingBox,pv,pv);
                %             elseif strcmp(type,'rectLong')
                %
                %                 pv=[0 0; 10 0; 10 1; 0 1; 0 0]/10;
                %                 [p,t]=distmesh2d(@dpoly,@huniform,el,boundingBox,pv,pv);
            elseif strcmp(type,'circle')
                fd=@(p) sqrt(sum(p.^2,2))-1/2;
                [p,t]=distmesh2d(fd,@huniform,el,[-1,-1;1,1],[]);
                p = p + 1/2;
            elseif strcmp(type,'rectangleCircularHole')
                fd=@(p) ddiff(drectangle(p,0,1,0,1),dcircle(p,0.5,0.5,0.25));
                fh=@(p) el+0.05*dcircle(p,0.5,0.5,0.25);
%                 fh=@(p) el;
                [p,t]=distmesh2d(fd,fh,el,boundingBox,[0,0;0,1;1,0;1,1]);
            elseif strcmp(type,'ellipse')
                fd=@(p) p(:,1).^2/2^2+p(:,2).^2/1^2-1;
                % el * 2 to account for scaling
                [p,t]=distmesh2d(fd,@huniform,el*2,[-2,-1;2,1],[]);
                p = p/2;
                p(:,1) = p(:,1) + 1;
                p(:,2) = p(:,2) + 0.5;
            elseif strcmp(type,'poly')
                
                pv=[-0.4 -0.5;0.4 -0.2;0.4 -0.7;1.5 -0.4;0.9 0.1;
                    1.6 0.8;0.5 0.5;0.2 1;0.1 0.4;-0.7 0.7;-0.4 -0.5];
                [p,t]=distmesh2d(@dpoly,@huniform,el*2,[-1,-1; 2,1],pv,pv);
                p = p/2;
                p(:,1) = p(:,1) - min(p(:,1));
                p(:,2) = p(:,2) - min(p(:,2));
                
            end
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
        trimesh(elem,nodeM);
        axis equal
        grid on
        % figure
        % simpplot(nodeM,elem)
        save([meshname '.mat'],'nodeM','elem')
        
    end
end

%% control the qualities of 2d tri mesh from distmesh using triangle with gptoolbox
% need to add external gptoolbox wrapper for triangle
% using mesh out put from distmesh as input to triangle
regenerate = false;

fs = filesep;

% type_list = {'triangle' 'circle' 'rectangleCircularHole' 'rect','poly','square' 'ellipse'};
type_list = {'triangle','circle','rect','rect_horizontal'};

% edge length list (parameters for distmesh)
% used to generate meshes from distmesh, which is then passed to
% triangle
el_list = [0.1 0.05];


% max area list (control parameter in triangle)
% maxA_list = [0.1 0.01 0.001 0.0005];
maxA_list = [0.01 0.001];
% now each area for quality control in triangle corresponds to each edge length from distmesh
% maxA is the actual knob to control the mesh


% everything in a fixed bounding box
boundingBox = [-0.1 -0.1; 1.1 1.1];

for i_type = 1:length(type_list)
    type = type_list{i_type};
    
    for i = 1:length(el_list)
        el = el_list(i);
        fprintf('%s\n',type);
        
        distmesh_meshname = sprintf('sim_data%ctri_meshes%c%s el%.e',fs,fs,type, el);
        %         meshname = sprintf('mesh_data%c%s_inv_edge_length_%.3d',fs,type, 1/el);
        
        if or(regenerate,exist([distmesh_meshname '.mat'], 'file') ~= 2)
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
            elseif strcmp(type,'rect_horizontal')
                pv=[0 0; 1 0; 1 3; 0 3; 0 0]/3;
                [p,t]=distmesh2d(@dpoly,@huniform,el,boundingBox,pv,pv);
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
            load([distmesh_meshname '.mat'], 'nodeM', 'elem');
        end
        save([distmesh_meshname '.mat'],'nodeM','elem')
        
    end
    
    for i = 1:length(maxA_list)
        el = el_list(i);
        maxA = maxA_list(i);
        fprintf('%s\n',type);
        distmesh_meshname = sprintf('sim_data%ctri_meshes%c%s el%.e',fs,fs,type, el);
        
        meshname = sprintf('mesh_data%c%s_maxA_%.d',fs,type, maxA);
        if or(regenerate,exist([meshname '.mat'], 'file') ~= 2)
            
            load([distmesh_meshname '.mat'], 'nodeM', 'elem');
            [TV,TF,TN] = triangle(nodeM, elem, 'MaxArea',maxA,'Quality',30);
            
            p = TV;
            t = TF;
            [p,t] = fixmesh(p,t);
            
            
            nodeM = p;
            elem = t;
            
        else
            load([meshname '.mat'], 'nodeM', 'elem');
        end
        colormap jet
        s_alpha = sprintf('%f',el);
        
        triplot(elem,nodeM(:,1),nodeM(:,2));
        axis equal
        grid on

        save([meshname '.mat'],'nodeM','elem')
        
        print(meshname,'-dpng')
    end
    
    
end

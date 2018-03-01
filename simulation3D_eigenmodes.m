
function out = simulation3D_eigenmodes(varargin)
close all
close(gcf)

out = 0;
% parse input

fs = filesep;

i_arg = 1;
while (i_arg <= nargin)
    switch varargin{i_arg}
        
        case 'dt'
            i_arg = i_arg + 1;
            dt = varargin{i_arg};
        case 'T'
            i_arg = i_arg + 1;
            T = varargin{i_arg};
        case 'meshfile'
            i_arg = i_arg + 1;
            meshfile = varargin{i_arg};
        case 'h'
            i_arg = i_arg + 1;
            h = varargin{i_arg};
            %         case 'simulation_type'
            %             i_arg = i_arg + 1;
            %             simulation_type = varargin{i_arg};
        case 'DGeta'
            i_arg = i_arg + 1;
            DGeta = varargin{i_arg};
        case 'solver'
            i_arg = i_arg + 1;
            solver = eval(['@' varargin{i_arg}]);
        case 'Y'
            i_arg = i_arg + 1;
            Y = varargin{i_arg};
        case 'P'
            i_arg = i_arg + 1;
            P = varargin{i_arg};
        case 'rho'
            i_arg = i_arg + 1;
            rho = varargin{i_arg};
        case 'a'
            i_arg = i_arg + 1;
            a = varargin{i_arg};
        case 'b'
            i_arg = i_arg + 1;
            b = varargin{i_arg};
        case 'material_type'
            i_arg = i_arg + 1;
            material_type = varargin{i_arg};
        case 'axis_box'
            i_arg = i_arg + 1;
            axis_box = varargin{i_arg};
        case 'constraint'
            i_arg = i_arg + 1;
            constraint = varargin{i_arg};
        case 'eig_modes'
            i_arg = i_arg + 1;
            eig_modes = varargin{i_arg};
        case 'gravity'
            i_arg = i_arg + 1;
            gravity = varargin{i_arg};
        case 'deformation_scale'
            i_arg = i_arg + 1;
            deformation_scale = varargin{i_arg};
        case 'reduced_dimension'
            i_arg = i_arg + 1;
            reduced_dimension = varargin{i_arg};
        case 'NMSC'
            i_arg = i_arg + 1;
            NMSC = varargin{i_arg};
        case 'fine_meshfile'
            i_arg = i_arg + 1;
            fine_meshfile = varargin{i_arg};
        case 'scene_name'
            i_arg = i_arg + 1;
            scene_name = varargin{i_arg};
            
    end
    i_arg = i_arg + 1;
end

tsteps = T/dt;
%%
scriptdir = script_directory(mfilename);
mkdir(scriptdir);
sim_directory_name = sim_directory(varargin);
simdir = [scriptdir,fs,sim_directory_name];
mkdir(simdir);

if (exist([simdir fs 'trajectory.mat'], 'file') ~= 2)
    
    if NMSC
        %% process the high res mesh first
        
        if exist([fine_meshfile '.mat'], 'file') ~= 2
            error('mesh does not exist')
        else
            load([fine_meshfile '.mat'], 'nodeM', 'elem');
            
        end
        %
        %                 if any(strfind(meshfile,'octopus'))
        %                     nodeM = nodeM(:,[3 2 1]);
        %                 end
        if any(strfind(meshfile,'armadillo'))
            nodeM(:,3) = -nodeM(:,3); % flip the armadillo
        end
        
        N = size(nodeM,1);
        %
        if any(strfind(meshfile,'small_bar'))
            [~, ind] = sortrows(nodeM,1);
            
            % put the nodes in order (in z)
            nodeM = nodeM(ind,:);
            
            % create the rank for the original z
            [~,R] = sort(ind);
            % the new element list
            elem = R(elem(:,:));
            
        elseif any(strfind(meshfile,'bunny'))
            [~, ind] = sortrows(nodeM,3);
            
            % put the nodes in order (in z)
            nodeM = nodeM(ind,:);
            
            % create the rank for the original z
            [~,R] = sort(ind);
            % the new element list
            elem = R(elem(:,:));
            
        elseif any(strfind(meshfile,'armadillo'))
            [~, ind] = sortrows(nodeM,3);
            
            % put the nodes in order (in z)
            nodeM = nodeM(ind,:);
            
            % create the rank for the original z
            [~,R] = sort(ind);
            % the new element list
            elem = R(elem(:,:));
            
        end
        
        % choose constraints and material for the scene
        [indLogical, Y_list, P_list, rho_list, constraint_points] = constraint_selection(meshfile, scene_name, nodeM, elem, Y, P, rho);
        
        obj = elasticTetObj(nodeM, elem);
        
        if strcmp(scene_name(1:4),'hete')
            switch material_type
                case 'linear'
                    obj.SetMaterial( Y_list, P_list, rho_list, 2 * ones(size(elem,1),1) , a, b); % set the tri to linear
                    %     obj.SetMaterial( Y, P, rho, 1:size(elem,1), 1); % set the tri to neo-hookean
                case 'neo-hookean'
                    obj.SetMaterial( Y_list, P_list, rho_list, 1 * ones(size(elem,1),1), a, b); % set the tri to neo-hookean
                case 'stvk'
                    obj.SetMaterial( Y_list, P_list, rho_list, 3 * ones(size(elem,1),1), a, b); % set the tri to stvk
            end
        else
            switch material_type
                case 'linear'
                    obj.SetMaterial( Y, P, rho, 2 , a, b); % set the tri to linear
                    %     obj.SetMaterial( Y, P, rho, 1:size(elem,1), 1); % set the tri to neo-hookean
                case 'neo-hookean'
                    obj.SetMaterial( Y, P, rho, 1, a, b); % set the tri to neo-hookean
                case 'stvk'
                    obj.SetMaterial( Y, P, rho, 3, a, b); % set the tri to stvk
            end
        end
        
        Dx = 0*rand(obj.Dim*N,1); % displacement field. set zero for rest state
        obj.SetCurrentState(Dx);
        %
        M = obj.M;
        K = obj.StiffnessMatrix;
        M = M(indLogical,indLogical);
        K = K(indLogical,indLogical);
        %         K = K(~ind_fix,~ind_fix); % extract the non-fixed part
        
        [V,D] = eigs(K,M,eig_modes,'sm');
        [low_eig, permutation_indices] = sort(diag(D));
        
        high_res_low_eig = low_eig;
        eigv = V(:,permutation_indices);
        
        fig = gcf;
        

        
        for i_mode = 1:eig_modes
%             Dx(indLogical) = eigv(:,i_mode)/20;
            Dx(indLogical) = eigv(:,i_mode);
            obj.SetCurrentState(Dx);
            delete(gca)
            ha = obj.init_vis;
            obj.simple_vis(obj.vis_handle);
            axis equal
            if any(strfind(meshfile,'small_bar'))
                campos([-1.1513   -1.6081    1.4854]);
                camtarget([0.0810   -0.0021    0.0680])
                camva(6.9295);
            elseif any(strfind(meshfile,'octopus'))
                campos([   -3.7639   -4.9106    3.4267]);
                camtarget([    0.0603    0.0732   -0.2002])
                camva(6.9295);
            elseif any(strfind(meshfile,'armadillo'))
                campos([2.4703  -20.8381    5.3614]);
                camtarget([0.7785   -0.6424    0.9054])
                camva(6.9295);
            elseif any(strfind(meshfile,'horse'))
%                 campos('auto')
%                 camtarget('auto')
                campos([    -12.6549   -9.8499   17.1693]);
                camtarget([    -0.0506   -0.0023    0.0166])
%                 camva(7.3921);
                camup([    0.5763    0.4503    0.6820]);
            end
            print(fig,[sim_directory_name(10:end) '-hi-res-mode' num2str(i_mode) '.png'],'-dpng')
            delete(gca)
            hold on
        end
    end
%     close(fig)
% cla
%%
    
    if exist([meshfile '.mat'], 'file') ~= 2
        error('mesh does not exist')
        
    else
        load([meshfile '.mat'], 'nodeM', 'elem');
        
    end
    
    %     if any(strfind(meshfile,'octopus'))
    %         nodeM = nodeM(:,[3 2 1]);
    %     end
    if any(strfind(meshfile,'armadillo'))
        nodeM(:,3) = -nodeM(:,3); % flip the armadillo
    end
    
    
    N = size(nodeM,1);
    
    
    
    if any(strfind(meshfile,'small_bar'))
        [~, ind] = sortrows(nodeM,1);
        
        % put the nodes in order (in z)
        nodeM = nodeM(ind,:);
        
        % create the rank for the original z
        [~,R] = sort(ind);
        % the new element list
        elem = R(elem(:,:));
        
    elseif any(strfind(meshfile,'bunny'))
        [~, ind] = sortrows(nodeM,3);
        
        % put the nodes in order (in z)
        nodeM = nodeM(ind,:);
        
        % create the rank for the original z
        [~,R] = sort(ind);
        % the new element list
        elem = R(elem(:,:));
        
    elseif any(strfind(meshfile,'armadillo'))
        [~, ind] = sortrows(nodeM,3);
        
        % put the nodes in order (in z)
        nodeM = nodeM(ind,:);
        
        % create the rank for the original z
        [~,R] = sort(ind);
        % the new element list
        elem = R(elem(:,:));
        
    end
    
    % choose constraints and material for the scene
    [indLogical, Y_list, P_list, rho_list, constraint_points] = constraint_selection(meshfile, scene_name, nodeM, elem, Y, P, rho);
    
    %             if any(strfind(meshfile,'armadillo'))
    %                 constraints = constraint_points;
    %             end
    
    obj = elasticTetObj(nodeM, elem);
    
    if strcmp(scene_name(1:4),'hete')
        switch material_type
            case 'linear'
                obj.SetMaterial( Y_list, P_list, rho_list, 2 * ones(size(elem,1),1) , a, b); % set the tri to linear
                %     obj.SetMaterial( Y, P, rho, 1:size(elem,1), 1); % set the tri to neo-hookean
            case 'neo-hookean'
                obj.SetMaterial( Y_list, P_list, rho_list, 1 * ones(size(elem,1),1), a, b); % set the tri to neo-hookean
            case 'stvk'
                obj.SetMaterial( Y_list, P_list, rho_list, 3 * ones(size(elem,1),1), a, b); % set the tri to stvk
        end
    else
        switch material_type
            case 'linear'
                obj.SetMaterial( Y, P, rho, 2, a, b); % set the tri to linear
                %     obj.SetMaterial( Y, P, rho, 1:size(elem,1), 1); % set the tri to neo-hookean
            case 'neo-hookean'
                obj.SetMaterial( Y, P, rho, 1, a, b); % set the tri to neo-hookean
            case 'stvk'
                obj.SetMaterial( Y, P, rho, 3, a, b); % set the tri to stvk
        end
    end
    
    switch gravity
        case 'on'
            obj.gravity_on = true;
            obj.calculateGravity;
        case 'off'
            obj.gravity_on = false;
            obj.calculateGravity;
    end
    
    Dx = 0*rand(obj.Dim*N,1);
    obj.SetCurrentState(Dx);
    M = obj.M;
    K = obj.StiffnessMatrix;
    M = M(indLogical,indLogical);
    K = K(indLogical,indLogical);
    %         K = K(~ind_fix,~ind_fix); % extract the non-fixed part
    
    [V,D] = eigs(K,M,eig_modes,'sm');
    [low_eig, permutation_indices] = sort(diag(D));
    
    obj.indLogical = indLogical;
    eigv = V(:,permutation_indices);
    

    v = zeros(length(Dx),1);
    u = [Dx; v];
    fig = gcf;

    % rate to draw the scene
    for i_mode = 1:eig_modes
%         Dx(indLogical) = eigv(:,i_mode)/20;
            Dx(indLogical) = eigv(:,i_mode);

        obj.SetCurrentState(Dx);
        delete(gca)
        ha = obj.init_vis;
        obj.simple_vis(obj.vis_handle);
        axis equal
        if any(strfind(meshfile,'small_bar'))
            campos([-1.1513   -1.6081    1.4854]);
            camtarget([0.0810   -0.0021    0.0680])
            camva(6.9295);
        elseif any(strfind(meshfile,'octopus'))
            campos([   -3.7639   -4.9106    3.4267]);
            camtarget([    0.0603    0.0732   -0.2002])
            camva(6.9295);
        elseif any(strfind(meshfile,'armadillo'))
            campos([2.4703  -20.8381    5.3614]);
            camtarget([0.7785   -0.6424    0.9054])
            camva(6.9295);
        elseif any(strfind(meshfile,'horse'))
                campos([    -12.6549   -9.8499   17.1693]);
                camtarget([    -0.0506   -0.0023    0.0166])
%                 camva(7.3921);
                camup([    0.5763    0.4503    0.6820]);
        end
        print(fig,[sim_directory_name(10:end) '-mode' num2str(i_mode) '.png'],'-dpng')
        delete(gca)
        hold on
    end
    
    
    
    
end
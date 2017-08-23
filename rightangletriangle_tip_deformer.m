function rightangletriangle_tip_deformer
% load mesh to analysis the eigenmodes
close all
clear all

draw = false;
rerun_flag = true;
% type_list = {'ellipse'}
% type_list = {'rectangleCircularHole'};
% type_list = {'rect'};

% max area list

maxA_list = [1 0.1 0.05 0.01]';

interpolation_step_list = [0.3];

fs = filesep;
mass_flag = true;

type = 'rightangletriangle';

axis_box = [-1 2 -1 2];
% options = optimoptions('fsolve','Algorithm','levenberg-marquardt');
options = optimoptions('fsolve','TolFun',1.e-9,'TolX',1.e-9,'Display','final');

energy_table = zeros(length(maxA_list),length(interpolation_step_list));

for i = 1:length(maxA_list)
    maxA = maxA_list(i);
    %         meshname = sprintf('simData%ctrimeshes%c%s maxA%.e',fs,fs,type, maxA);
    
    meshname = sprintf('simData%ctrimeshes%c%s maxA%.e tritool%s',fs,fs,type, maxA);
    if exist([meshname '.mat'], 'file') ~= 2
        disp('mesh does not exist')
        break
    else
        load([meshname '.mat'], 'nodeM', 'elem');
        
    end
    
    
    elem(:,[1 3]) = elem(:,[3 1]);
    
    N = size(nodeM,1);
    
    
    
    Y = 5; % Young's modululs
    P = 0.4; % Poisson ratio
    rho = 1; % density
    
    filename = sprintf('simData%ctrimeshes%cmatlabData%c%s maxA%.e mass%d Y%.e P%.e rho%.e', fs, fs, fs, type, maxA, mass_flag, Y, P, rho);
    if (exist([filename '.mat'], 'file') ~= 2) | rerun_flag
        
        % construct triangular mesh object
        obj = elasticTriObj(nodeM, elem);
        
        
        obj.SetMaterial( Y, P, rho, 1:size(elem,1), 1); % set the tri to neohookean
        
%         M = obj.M;
        % print out the total mass for sanity check
        % should be the same for all mesh resolution
        % disp('total mass')
        % disp(sum(spdiags(M)))
        
%         M = M(~ind_fix,~ind_fix); % extract the non-fixed part
%         
        Dx = 0*rand(2*N,1); % displacement field. set zero for rest state
        obj.SetCurrentState(Dx);
%         
%         K = obj.StiffnessMatrix;
%         K = K(~ind_fix,~ind_fix); % extract the non-fixed part
        
        save([filename '.mat']);
    else
        load([filename '.mat'], 'obj');
    end
    
    
    for i_deform = 1:length(interpolation_step_list)
        interpolation_step = interpolation_step_list(i_deform);
        if maxA == 1
            [~,ind_apex] = max(nodeM(:,2));
            Dx = zeros(size(obj.x));
            Dx(2*(ind_apex-1)+1 : 2*ind_apex) = [-1 -1]';
            defo = Dx * interpolation_step;
            obj.SetCurrentState(defo)
            node = obj.node;
            trimesh(elem,node(:,1),node(:,2));
            
            hold on
            ef = obj.ElasticForce;

            quiver(node(:,1),node(:,2),ef(1:2:end),ef(2:2:end),0)
            hold off

            
        else
            Xind_y_axis = abs(nodeM(:,1)) < 1e-3;
            ind_y_axis = reshape(transpose(repmat(Xind_y_axis,1,2)),[],1); % logical index for total position vector
            ratio_y_axis = nodeM(Xind_y_axis,2);
            
            Dx_final_y_axis = [((-1 - 0) * ratio_y_axis -0), (0 - nodeM(Xind_y_axis,2))];
            Dx_final_y_axis = reshape(transpose(Dx_final_y_axis),[],1);
            
            Xind_slope = (abs(nodeM(:,1)+nodeM(:,2) - 1 ) < 1e-3) & abs(nodeM(:,1) < 1);
            ind_slope = reshape(transpose(repmat(Xind_slope,1,2)),[],1); % logical index for total position vector
            ratio_slope = nodeM(Xind_slope,2);
            Dx_final_slope = [((-1 -1) * ratio_slope + 1) - nodeM(Xind_slope,1), (0 - nodeM(Xind_slope,2))];
            Dx_final_slope = reshape(transpose(Dx_final_slope),[],1);
            
            
            
            Xind_boundary = (abs(nodeM(:,1)) < 1e-3) | (abs(nodeM(:,1)+nodeM(:,2) - 1 ) < 1e-3);
            ind_boundary = reshape(transpose(repmat(Xind_boundary,1,2)),[],1); % logical index for total position vector
            
            Xind_apex = nodeM(:,2) == max(nodeM(:,2));
            ind_apex = reshape(transpose(repmat(Xind_apex,1,2)),[],1); % logical index for total position vector
            Xind_apex2 = (nodeM(:,1) == max(nodeM(:,1))) & (nodeM(:,2) == min(nodeM(:,2)));
            ind_apex2 = reshape(transpose(repmat(Xind_apex2,1,2)),[],1); % logical index for total position vector
            Xind_apex3 = (nodeM(:,1) == min(nodeM(:,1))) & (nodeM(:,2) == min(nodeM(:,2)));
            ind_apex3 = reshape(transpose(repmat(Xind_apex3,1,2)),[],1); % logical index for total position vector
            
            
            
%             Xind_fix = abs(nodeM(:,2)) < 1e-3;
%             ind_fix = reshape(transpose(repmat(Xind_fix,1,2)),[],1); % logical index for total position vector
            ind_fix = ind_apex | ind_apex2 | ind_apex3;
%             Dx = Dx_warm_start;
            Dx = zeros(size(obj.X));
            Dx(ind_apex) = [-1 -1]';
            Dx(ind_y_axis) = Dx_final_y_axis;
            Dx(ind_slope) = Dx_final_slope;
            
            defo = Dx * interpolation_step;
            obj.SetCurrentState(defo);
            
            P = [0 0;...
                0 1;...
                1 0];
            T = [1 2 3];
            TR = triangulation(T,P);
            defo_warm_start = zeros(size(defo));
            for i_n = 1:obj.N
                PC = obj.X(2*(i_n-1)+1:2*i_n)';
                B = cartesianToBarycentric(TR,1,PC);
                defo_warm_start(2*(i_n-1)+1:2*i_n) = B(1) * [0 0]' + B(2) * obj.x(ind_apex) + B(3) * [1 0]' - PC';
            end
            
            free_defo = defo_warm_start(~ind_fix);
            dx = fsolve(@nonlinear_f, free_defo, options);
            
            defo(~ind_fix) = dx; 
            obj.SetCurrentState(defo)
            node = obj.node;
            figure
            trimesh(elem,node(:,1),node(:,2));
            hold on
            ef = obj.ElasticForce;

            quiver(node(:,1),node(:,2),ef(1:2:end),ef(2:2:end),0)
            hold off
        end
        
        axis(axis_box)
        drawnow
        print([filename num2str(i_deform)],'-dpng');
        energy_table(i, i_deform) = obj.totalEnergy;
    end
    
end

    array2table(energy_table, 'RowNames', matlab.lang.makeValidName(cellstr((num2str(maxA_list)))),...
        'VariableNames', matlab.lang.makeValidName(cellstr(num2str(interpolation_step_list'))))
    
    function out = nonlinear_f(free_dx)
        defo(~ind_fix) = free_dx; 
        obj.SetCurrentState(defo);
        force =  obj.ElasticForce;
        out = force(~ind_fix);
    end

end

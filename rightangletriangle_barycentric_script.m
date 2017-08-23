function rightangletriangle_barycentric_script
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


fs = filesep;

type = 'rightangletriangle';




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
    
    filename = sprintf('simData%ctrimeshes%cmatlabData%c%s maxA%.e barycentric', fs, fs, fs, type, maxA);
    
    
    Xind_apex = nodeM(:,2) == max(nodeM(:,2));
    ind_apex = reshape(transpose(repmat(Xind_apex,1,2)),[],1); % logical index for total position vector
    Xind_apex2 = (nodeM(:,1) == max(nodeM(:,1))) & (nodeM(:,2) == min(nodeM(:,2)));
    ind_apex2 = reshape(transpose(repmat(Xind_apex2,1,2)),[],1); % logical index for total position vector
    Xind_apex3 = (nodeM(:,1) == min(nodeM(:,1))) & (nodeM(:,2) == min(nodeM(:,2)));
    ind_apex3 = reshape(transpose(repmat(Xind_apex3,1,2)),[],1); % logical index for total position vector
    
    Xind_fix = Xind_apex | Xind_apex2 | Xind_apex3;
    ind_fix = ind_apex | ind_apex2 | ind_apex3;
    T = 1:3;
    
    TR = triangulation(T,nodeM(Xind_fix,1),nodeM(Xind_fix,2));
    Bary = cartesianToBarycentric(TR,ones(N,1),nodeM);
    
    save([filename '.mat'], 'Bary');
end
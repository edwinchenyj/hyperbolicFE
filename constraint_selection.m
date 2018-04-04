function [indLogical, Y_list, P_list, rho_list, constraint_points] = constraint_selection(meshfile, scene_name, nodeM, elem, Y, P, rho)
Y_list = Y * ones(size(elem,1),1);
P_list = P * ones(size(elem,1),1);
rho_list = rho * ones(size(elem,1),1);
constraint_points = 0;
if any(strfind(meshfile,'small_bar'))
    switch scene_name
        case 'hete1'
            % constraint more points on the left
            left_points = find(abs(nodeM(:,1)-0) < 0.1/3);
            right_points = find(abs(nodeM(:,1)-1/5) < 0.01/5);
            
            % get logical indices for the left and right nodes
            
            indLeft = left_points;
            indRight = right_points;
            
            % nFixed = length(indRight);
            
            N = size(nodeM,1); % number of nodes
            indAll = 1:N;
            indRemove = indAll([indLeft(:)]);
            indLogical = logical(ones(3,N));
            indLogical(:,indRemove) = logical(0);
            indLogical = indLogical(:);
            
            % %
            % % find elements close to constraints and make them softer
            soft_points = find(abs(nodeM(:,1)-0) > 1/5/2);
            
            soft_points_height_threshold = min(soft_points);
            
            Y_list = Y * ones(size(elem,1),1);
            P_list = P * ones(size(elem,1),1);
            rho_list = rho * ones(size(elem,1),1);
            %
            for i = 1:size(elem,1)
                if all(elem(i,:) > soft_points_height_threshold)
                    Y_list(i) = Y_list(i) * 10;
                end
            end
        case 'unif1'
            % constraint more points on the left
            left_points = find(abs(nodeM(:,1)-0) < 0.01/3);
            right_points = find(abs(nodeM(:,1)-1/5) < 0.01/5);
            
            % get logical indices for the left and right nodes
            
            indLeft = left_points;
            indRight = right_points;
            
            % nFixed = length(indRight);
            
            N = size(nodeM,1); % number of nodes
            indAll = 1:N;
            indRemove = indAll([indLeft(:)]);
            indLogical = logical(ones(3,N));
            indLogical(:,indRemove) = logical(0);
            indLogical = indLogical(:);
    end
elseif any(strfind(meshfile,'barmesh'))
    switch scene_name
        case 'hete1'
            % constraint more points on the left
            left_points = find(abs(nodeM(:,1)-min(nodeM(:,1))) < 0.01/3);
            right_points = find(abs(nodeM(:,1)-max(nodeM(:,1))) < 0.01/5);
            
            % get logical indices for the left and right nodes
            
            indLeft = left_points;
            indRight = right_points;
            
            % nFixed = length(indRight);
            
            N = size(nodeM,1); % number of nodes
            indAll = 1:N;
            indRemove = indAll([indRight(:)]);
            indLogical = logical(ones(3,N));
            indLogical(:,indRemove) = logical(0);
            indLogical = indLogical(:);
            
            constraint_points = right_points;
            
            % %
            % % find elements close to constraints and make them softer
            soft_points = find(abs(nodeM(:,1)-0) > 1/5*2);
            
            soft_points_height_threshold = min(soft_points);
            
            Y_list = Y * ones(size(elem,1),1);
            P_list = P * ones(size(elem,1),1);
            rho_list = rho * ones(size(elem,1),1);
            %
            for i = 1:size(elem,1)
                if all(elem(i,:) > soft_points_height_threshold)
                    Y_list(i) = Y_list(i) * 10;
                end
            end
        case 'hete2'
            % constraint more points on the left
            left_points = find(abs(nodeM(:,1)-min(nodeM(:,1))) < 0.01/3);
            right_points = find(abs(nodeM(:,1)-max(nodeM(:,1))) < 0.01/5);
            
            % get logical indices for the left and right nodes
            
            indLeft = left_points;
            indRight = right_points;
            
            % nFixed = length(indRight);
            
            N = size(nodeM,1); % number of nodes
            indAll = 1:N;
            indRemove = indAll([indRight(:)]);
            indLogical = logical(ones(3,N));
            indLogical(:,indRemove) = logical(0);
            indLogical = indLogical(:);
            
            constraint_points = right_points;
            
            % %
            % % find elements close to constraints and make them softer
            soft_points = find(abs(nodeM(:,1)-max(nodeM(:,1))) < 1/5*2);
            
            soft_points_height_threshold = min(soft_points);
            
            Y_list = Y * ones(size(elem,1),1);
            P_list = P * ones(size(elem,1),1);
            rho_list = rho * ones(size(elem,1),1);
            %
            for i = 1:size(elem,1)
                if all(elem(i,:) > soft_points_height_threshold)
                    Y_list(i) = Y_list(i)/2;
                end
            end
        case 'unif1'
            
            % constraint more points on the left
            left_points = find(abs(nodeM(:,1)-min(nodeM(:,1))) < 0.01/3);
            right_points = find(abs(nodeM(:,1)-max(nodeM(:,1))) < 0.01/5);
            
            % get logical indices for the left and right nodes
            
            indLeft = left_points;
            indRight = right_points;
            
            % nFixed = length(indRight);
            
            N = size(nodeM,1); % number of nodes
            indAll = 1:N;
            indRemove = indAll([indRight(:)]);
            indLogical = logical(ones(3,N));
            indLogical(:,indRemove) = logical(0);
            indLogical = indLogical(:);
            
            constraint_points = right_points;
    end
elseif any(strfind(meshfile,'octopus'))
    switch scene_name
        case 'hete1'
        case 'unif1'
            center_points = find(and(and(abs(nodeM(:,1)-0) < 0.1, ...
                abs(nodeM(:,2)+0.1) < 0.1),abs(nodeM(:,3)- 0.1) < 0.1));
            indCenter = center_points;
            N = size(nodeM,1); % number of nodes
            indAll = 1:N;
            indRemove = indAll([indCenter(:)]);
            indLogical = logical(ones(3,N));
            indLogical(:,indRemove) = logical(0);
            indLogical = indLogical(:);
            
        case 'unif2'
            center_points = find(and(and(abs(nodeM(:,1)-0) < 0.15, ...
                abs(nodeM(:,2)+0.1) < 0.15),abs(nodeM(:,3)- 0.1) < 0.15));
            indCenter = center_points;
            N = size(nodeM,1); % number of nodes
            indAll = 1:N;
            indRemove = indAll([indCenter(:)]);
            indLogical = logical(ones(3,N));
            indLogical(:,indRemove) = logical(0);
            indLogical = indLogical(:);
            
        case 'unif3'
            center_points = find(and(and(abs(nodeM(:,1)-0) < 0.2, ...
                abs(nodeM(:,2)+0.1) < 0.2),abs(nodeM(:,3)- 0.1) < 0.2));
            indCenter = center_points;
            N = size(nodeM,1); % number of nodes
            indAll = 1:N;
            indRemove = indAll([indCenter(:)]);
            indLogical = logical(ones(3,N));
            indLogical(:,indRemove) = logical(0);
            indLogical = indLogical(:);
        case 'unif4'
            center_points = find(and(and(abs(nodeM(:,1)-0.22) < 1.1, ...
                abs(nodeM(:,2)+1.1865) < 1.1135),abs(nodeM(:,3)- 1.1) < 2.2));
            indCenter = center_points;
            N = size(nodeM,1); % number of nodes
            indAll = 1:N;
            indRemove = indAll([indCenter(:)]);
            indLogical = logical(ones(3,N));
            indLogical(:,indRemove) = logical(0);
            indLogical = indLogical(:);
            
    end
elseif any(strfind(meshfile,'armadillo'))
    switch scene_name
        case 'hete1'
            bottom_points = find(abs(nodeM(:,3)-min(nodeM(:,3))) < 0.08);
            top_points = find(abs(nodeM(:,3)-max(nodeM(:,3))) < 0.08);
            
            %% get logical indices
            indBottom = bottom_points;
            indTop = top_points;
            constraint_points = top_points;
            nFixed = length(indTop);
            
            N = size(nodeM,1); % number of nodes
            indAll = 1:N;
            indRemove = indAll([indTop(:)]);
            indLogical = logical(ones(3,N));
            indLogical(:,indRemove) = logical(0);
            indLogical = indLogical(:);
            
            
            %% find elements close to constraints and make them softer
            soft_points = find(abs(nodeM(:,3)-max(nodeM(:,3))) < 0.2);
            
            soft_points_height_threshold = min(soft_points);
            
            Y_list = Y * ones(size(elem,1),1);
            P_list = P * ones(size(elem,1),1);
            rho_list = rho * ones(size(elem,1),1);
            
            for i = 1:size(elem,1)
                if all(elem(i,:) > soft_points_height_threshold)
                    Y_list(i) = Y_list(i) / 10;
                end
            end
        case 'unif1'
            bottom_points = find(abs(nodeM(:,3)-min(nodeM(:,3))) < 0.08);
            top_points = find(abs(nodeM(:,3)-max(nodeM(:,3))) < 0.08);
            
            %% get logical indices
            indBottom = bottom_points;
            indTop = top_points;
            constraint_points = top_points;
            nFixed = length(indTop);
            
            N = size(nodeM,1); % number of nodes
            indAll = 1:N;
            indRemove = indAll([indTop(:)]);
            indLogical = logical(ones(3,N));
            indLogical(:,indRemove) = logical(0);
            indLogical = indLogical(:);
            
    end
    
elseif any(strfind(meshfile,'bunnymesh'))
    switch scene_name
        case 'hete1'
            bottom_points = find(abs(nodeM(:,3)-min(nodeM(:,3))) < 0.08);
            top_points = find(abs(nodeM(:,3)-max(nodeM(:,3))) < 0.08);
            
            %% get logical indices
            indBottom = bottom_points;
            indTop = top_points;
            constraint_points = top_points;
            nFixed = length(indTop);
            
            N = size(nodeM,1); % number of nodes
            indAll = 1:N;
            indRemove = indAll([indTop(:)]);
            indLogical = logical(ones(3,N));
            indLogical(:,indRemove) = logical(0);
            indLogical = indLogical(:);
            
            
            %% find elements close to constraints and make them softer
            soft_points = find(abs(nodeM(:,3)-max(nodeM(:,3))) < 0.2);
            
            soft_points_height_threshold = min(soft_points);
            
            Y_list = Y * ones(size(elem,1),1);
            P_list = P * ones(size(elem,1),1);
            rho_list = rho * ones(size(elem,1),1);
            
            for i = 1:size(elem,1)
                if all(elem(i,:) > soft_points_height_threshold)
                    Y_list(i) = Y_list(i) / 10;
                end
            end
        case 'hete2'
            bottom_points = find(abs(nodeM(:,3)-min(nodeM(:,3))) < 0.3);
            top_points = find(abs(nodeM(:,3)-max(nodeM(:,3))) < 0.3);
            
            %% get logical indices
            indBottom = bottom_points;
            indTop = top_points;
            constraint_points = top_points;
            nFixed = length(indTop);
            
            N = size(nodeM,1); % number of nodes
            indAll = 1:N;
            indRemove = indAll([indTop(:)]);
            indLogical = logical(ones(3,N));
            indLogical(:,indRemove) = logical(0);
            indLogical = indLogical(:);
            
            
            %% find elements close to constraints and make them softer
            soft_points = find(abs(nodeM(:,3)-max(nodeM(:,3))) < 0.5);
            
            soft_points_height_threshold = min(soft_points);
            
            Y_list = Y * ones(size(elem,1),1);
            P_list = P * ones(size(elem,1),1);
            rho_list = rho * ones(size(elem,1),1);
            
            for i = 1:size(elem,1)
                if all(elem(i,:) > soft_points_height_threshold)
                    Y_list(i) = Y_list(i) / 10;
                end
            end
        case 'hete3'
            bottom_points = find(abs(nodeM(:,3)-min(nodeM(:,3))) < 0.08);
            top_points = find(abs(nodeM(:,3)-max(nodeM(:,3))) < 0.08);
            
            %% get logical indices
            indBottom = bottom_points;
            indTop = top_points;
            constraint_points = top_points;
            nFixed = length(indTop);
            
            N = size(nodeM,1); % number of nodes
            indAll = 1:N;
            indRemove = indAll([indTop(:)]);
            indLogical = logical(ones(3,N));
            indLogical(:,indRemove) = logical(0);
            indLogical = indLogical(:);
            
            
            %% find elements close to constraints and make them softer
            soft_points = find(abs(nodeM(:,3)-max(nodeM(:,3))) < 0.2);
            
            soft_points_height_threshold = min(soft_points);
            
            Y_list = Y * ones(size(elem,1),1);
            P_list = P * ones(size(elem,1),1);
            rho_list = rho * ones(size(elem,1),1);
            
            for i = 1:size(elem,1)
                if all(elem(i,:) > soft_points_height_threshold)
                    Y_list(i) = Y_list(i) / 10;
                end
            end
        case 'hete4'
            bottom_points = find(abs(nodeM(:,3)-min(nodeM(:,3))) < 0.15);
            top_points = find(abs(nodeM(:,3)-max(nodeM(:,3))) < 0.15);
            
            %% get logical indices
            indBottom = bottom_points;
            indTop = top_points;
            constraint_points = top_points;
            nFixed = length(indTop);
            
            N = size(nodeM,1); % number of nodes
            indAll = 1:N;
            indRemove = indAll([indTop(:)]);
            indLogical = logical(ones(3,N));
            indLogical(:,indRemove) = logical(0);
            indLogical = indLogical(:);
            
            
            %% find elements close to constraints and make them softer
            soft_points = find(abs(nodeM(:,3)-max(nodeM(:,3))) < 0.2);
            
            soft_points_height_threshold = min(soft_points);
            
            Y_list = Y * ones(size(elem,1),1);
            P_list = P * ones(size(elem,1),1);
            rho_list = rho * ones(size(elem,1),1);
            
            for i = 1:size(elem,1)
                if all(elem(i,:) > soft_points_height_threshold)
                    Y_list(i) = Y_list(i) / 10;
                end
            end
        case 'hete5'
            
% fidn all the points around the same height with the lowest point
bottom_points = find(abs(nodeM(:,3)-nodeM(1,3)) < (0.01));

nFixed = length(bottom_points);

N = size(nodeM,1); % number of nodes

% N = size(nodeM,1); % number of nodes
indAll = 1:N;
indRemove = indAll(bottom_points(:));
indLogical = true(3,N);
indLogical(:,indRemove) = false;
indLogical = indLogical(:);
%% find elements close to constraints and make them softer
soft_points = find(abs(nodeM(:,3)-min(nodeM(1,3))) < (0.4));

soft_points_height_threshold = max(soft_points);

Y_list = Y * ones(size(elem,1),1);
P_list = P * ones(size(elem,1),1);
rho_list = rho * ones(size(elem,1),1);

for i = 1:size(elem,1)
    if all(elem(i,:) < soft_points_height_threshold)
        Y_list(i) = Y_list(i)/10;
    end
end
        case 'unif1'
            % constraint more points on the left
            top_points = find(abs(nodeM(:,3)-max(nodeM(:,3))) < 0.12);
            
            %% get logical indices
            %             indBottom = bottom_points;
            indTop = top_points;
            constraint_points = top_points;
            nFixed = length(indTop);
            
            N = size(nodeM,1); % number of nodes
            indAll = 1:N;
            indRemove = indAll([indTop(:)]);
            indLogical = logical(ones(3,N));
            indLogical(:,indRemove) = logical(0);
            indLogical = indLogical(:);
        case 'unif2'
            % constraint more points on the left
            top_points = find(abs(nodeM(:,3)-max(nodeM(:,3))) < 0.08);
            
            %% get logical indices
            %             indBottom = bottom_points;
            indTop = top_points;
            constraint_points = top_points;
            nFixed = length(indTop);
            
            N = size(nodeM,1); % number of nodes
            indAll = 1:N;
            indRemove = indAll([indTop(:)]);
            indLogical = logical(ones(3,N));
            indLogical(:,indRemove) = logical(0);
            indLogical = indLogical(:);
        case 'unif3'
            % constraint more points on the left
            top_points = find(abs(nodeM(:,3)-max(nodeM(:,3))) < 0.3);
            
            %% get logical indices
            %             indBottom = bottom_points;
            indTop = top_points;
            constraint_points = top_points;
            nFixed = length(indTop);
            
            N = size(nodeM,1); % number of nodes
            indAll = 1:N;
            indRemove = indAll([indTop(:)]);
            indLogical = logical(ones(3,N));
            indLogical(:,indRemove) = logical(0);
            indLogical = indLogical(:);
        case 'unif4'
            
% fidn all the points around the same height with the lowest point
bottom_points = find(abs(nodeM(:,3)-nodeM(1,3)) < (0.01));

nFixed = length(bottom_points);

N = size(nodeM,1); % number of nodes

% N = size(nodeM,1); % number of nodes
indAll = 1:N;
indRemove = indAll(bottom_points(:));
indLogical = true(3,N);
indLogical(:,indRemove) = false;
indLogical = indLogical(:);
    end
elseif any(strfind(meshfile,'horse'))
    switch scene_name
        %         case 'hete1'
        case 'unif1'
            % constraint more points on the left
            left_points = find(nodeM(:,2) > 0.06);
            
            % get logical indices for the left and right nodes
            
            indLeft = left_points;
            %             indRight = right_points;
            
            % nFixed = length(indRight);
            
            N = size(nodeM,1); % number of nodes
            indAll = 1:N;
            indRemove = indAll([indLeft(:)]);
            indLogical = logical(ones(3,N));
            indLogical(:,indRemove) = logical(0);
            indLogical = indLogical(:);
            
            constraint_points = left_points;
        case 'unif2'
            % constraint more points on the left
            Zmin = min(nodeM(:,3));
            bottom_points = find(nodeM(:,3) - Zmin < 0.05);
            
            % get logical indices for the left and right nodes
            
            indBottom = bottom_points;
            %             indRight = right_points;
            
            % nFixed = length(indRight);
            
            N = size(nodeM,1); % number of nodes
            indAll = 1:N;
            indRemove = indAll([indBottom(:)]);
            indLogical = logical(ones(3,N));
            indLogical(:,indRemove) = logical(0);
            indLogical = indLogical(:);
            
            constraint_points = bottom_points;
        case 'unif3'
            % constraint more points on the left
            Zmax = max(nodeM(:,3));
            top_points = find(abs(nodeM(:,3) - Zmax) < 0.036);
            
            % get logical indices for the left and right nodes
            
            indTop = top_points;
            %             indRight = right_points;
            
            % nFixed = length(indRight);
            
            N = size(nodeM,1); % number of nodes
            indAll = 1:N;
            indRemove = indAll([indTop(:)]);
            indLogical = logical(ones(3,N));
            indLogical(:,indRemove) = logical(0);
            indLogical = indLogical(:);
            
            constraint_points = top_points;
        case 'unif4'
            % constraint more points on the left
            Ymax = max(nodeM(:,2));
            right_points = find(abs(nodeM(:,2) - Ymax) < 0.1);
            
            % get logical indices for the left and right nodes
            
            %             right_points;
            %             indRight = right_points;
            
            % nFixed = length(indRight);
            
            N = size(nodeM,1); % number of nodes
            indAll = 1:N;
            indRemove = indAll([right_points(:)]);
            indLogical = logical(ones(3,N));
            indLogical(:,indRemove) = logical(0);
            indLogical = indLogical(:);
            
            constraint_points = right_points;
        case 'unif5'
            % constraint more points on the left
            top_points = find(abs(nodeM(:,3)-max(nodeM(:,3))) < 0.5/100);
            
            %% get logical indices
            %             indBottom = bottom_points;
            indTop = top_points;
            constraint_points = top_points;
            nFixed = length(indTop);
            
            N = size(nodeM,1); % number of nodes
            indAll = 1:N;
            indRemove = indAll([indTop(:)]);
            indLogical = logical(ones(3,N));
            indLogical(:,indRemove) = logical(0);
            indLogical = indLogical(:);
        case 'unif6'
            % constraint more points on the left
            top_points = find(abs(nodeM(:,3)-max(nodeM(:,3))) < 0.12);
            
            %% get logical indices
            %             indBottom = bottom_points;
            indTop = top_points;
            constraint_points = top_points;
            nFixed = length(indTop);
            
            N = size(nodeM,1); % number of nodes
            indAll = 1:N;
            indRemove = indAll([indTop(:)]);
            indLogical = logical(ones(3,N));
            indLogical(:,indRemove) = logical(0);
            indLogical = indLogical(:);
            
        case 'hete1'
            top_points = find(abs(nodeM(:,3)-max(nodeM(:,3))) < 0.2/100);
            
            %% get logical indices
            %             indBottom = bottom_points;
            indTop = top_points;
            constraint_points = top_points;
            nFixed = length(indTop);
            
            N = size(nodeM,1); % number of nodes
            indAll = 1:N;
            indRemove = indAll([indTop(:)]);
            indLogical = logical(ones(3,N));
            indLogical(:,indRemove) = logical(0);
            indLogical = indLogical(:);
            
            
            %% find elements close to constraints and make them softer
            soft_points = find(abs(nodeM(:,3)-max(nodeM(:,3))) < 0.3/100);
            
            soft_points_height_threshold = min(soft_points);
            
            Y_list = Y * ones(size(elem,1),1);
            P_list = P * ones(size(elem,1),1);
            rho_list = rho * ones(size(elem,1),1);
            
            for i = 1:size(elem,1)
                if all(elem(i,:) > soft_points_height_threshold)
                    Y_list(i) = Y_list(i) / 10;
                end
            end
        case 'hete2'
            top_points = find(abs(nodeM(:,3)-max(nodeM(:,3))) < 0.5/100);
            
            %% get logical indices
            %             indBottom = bottom_points;
            indTop = top_points;
            constraint_points = top_points;
            nFixed = length(indTop);
            
            N = size(nodeM,1); % number of nodes
            indAll = 1:N;
            indRemove = indAll([indTop(:)]);
            indLogical = logical(ones(3,N));
            indLogical(:,indRemove) = logical(0);
            indLogical = indLogical(:);
            
            
            %% find elements close to constraints and make them softer
            soft_points = find(abs(nodeM(:,3)-max(nodeM(:,3))) < 1/100);
            
            soft_points_height_threshold = min(soft_points);
            
            Y_list = Y * ones(size(elem,1),1);
            P_list = P * ones(size(elem,1),1);
            rho_list = rho * ones(size(elem,1),1);
            
            for i = 1:size(elem,1)
                if all(elem(i,:) > soft_points_height_threshold)
                    Y_list(i) = Y_list(i) / 10;
                end
            end
        case 'hete3'
            top_points = find(abs(nodeM(:,3)-max(nodeM(:,3))) < 0.12);
            
            %% get logical indices
            %             indBottom = bottom_points;
            indTop = top_points;
            constraint_points = top_points;
            nFixed = length(indTop);
            
            N = size(nodeM,1); % number of nodes
            indAll = 1:N;
            indRemove = indAll([indTop(:)]);
            indLogical = logical(ones(3,N));
            indLogical(:,indRemove) = logical(0);
            indLogical = indLogical(:);
            
            
            %% find elements close to constraints and make them softer
            soft_points = find(abs(nodeM(:,3)-max(nodeM(:,3))) < 0.25);
            
            soft_points_height_threshold = min(soft_points);
            
            Y_list = Y * ones(size(elem,1),1);
            P_list = P * ones(size(elem,1),1);
            rho_list = rho * ones(size(elem,1),1);
            
            for i = 1:size(elem,1)
                if all(elem(i,:) > soft_points_height_threshold)
                    Y_list(i) = Y_list(i) / 10;
                end
            end
        case 'hete4'
            top_points = find(abs(nodeM(:,3)-max(nodeM(:,3))) < 0.4);
            
            %% get logical indices
            %             indBottom = bottom_points;
            indTop = top_points;
            constraint_points = top_points;
            nFixed = length(indTop);
            
            N = size(nodeM,1); % number of nodes
            indAll = 1:N;
            indRemove = indAll([indTop(:)]);
            indLogical = logical(ones(3,N));
            indLogical(:,indRemove) = logical(0);
            indLogical = indLogical(:);
            
            
            %% find elements close to constraints and make them softer
            soft_points = find(abs(nodeM(:,3)-max(nodeM(:,3))) < 0.6);
            
            soft_points_height_threshold = min(soft_points);
            
            Y_list = Y * ones(size(elem,1),1);
            P_list = P * ones(size(elem,1),1);
            rho_list = rho * ones(size(elem,1),1);
            
            for i = 1:size(elem,1)
                if all(elem(i,:) > soft_points_height_threshold)
                    Y_list(i) = Y_list(i) / 10;
                end
            end
        case 'hete5'
            top_points = find(abs(nodeM(:,3)-max(nodeM(:,3))) < 0.1);
            
            %% get logical indices
            %             indBottom = bottom_points;
            indTop = top_points;
            constraint_points = top_points;
            nFixed = length(indTop);
            
            N = size(nodeM,1); % number of nodes
            indAll = 1:N;
            indRemove = indAll([indTop(:)]);
            indLogical = logical(ones(3,N));
            indLogical(:,indRemove) = logical(0);
            indLogical = indLogical(:);
            
            
            %% find elements close to constraints and make them softer
            soft_points = find(abs(nodeM(:,3)-max(nodeM(:,3))) < 0.3);
            
            soft_points_height_threshold = min(soft_points);
            
            Y_list = Y * ones(size(elem,1),1);
            P_list = P * ones(size(elem,1),1);
            rho_list = rho * ones(size(elem,1),1);
            
            for i = 1:size(elem,1)
                if all(elem(i,:) > soft_points_height_threshold)
                    Y_list(i) = Y_list(i) / 10;
                end
            end
        case 'hete6'
            top_points = find(abs(nodeM(:,3)-max(nodeM(:,3))) < 0.1);
            
            %% get logical indices
            %             indBottom = bottom_points;
            indTop = top_points;
            constraint_points = top_points;
            nFixed = length(indTop);
            
            N = size(nodeM,1); % number of nodes
            indAll = 1:N;
            indRemove = indAll([indTop(:)]);
            indLogical = logical(ones(3,N));
            indLogical(:,indRemove) = logical(0);
            indLogical = indLogical(:);
            
            
            %% find elements close to constraints and make them softer
            soft_points = find(abs(nodeM(:,3)-max(nodeM(:,3))) < 0.25);
            
            soft_points_height_threshold = min(soft_points);
            
            Y_list = Y * ones(size(elem,1),1);
            P_list = P * ones(size(elem,1),1);
            rho_list = rho * ones(size(elem,1),1);
            
            for i = 1:size(elem,1)
                if all(elem(i,:) > soft_points_height_threshold)
                    Y_list(i) = Y_list(i) / 5;
                end
            end
    end
    
elseif any(strfind(meshfile,'torusmesh'))
    switch scene_name
        case 'hete1'
            bottom_points = find(abs(nodeM(:,3)-min(nodeM(:,3))) < 0.08);
            top_points = find(abs(nodeM(:,3)-max(nodeM(:,3))) < 0.08);
            
            %% get logical indices
            indBottom = bottom_points;
            indTop = top_points;
            constraint_points = top_points;
            nFixed = length(indTop);
            
            N = size(nodeM,1); % number of nodes
            indAll = 1:N;
            indRemove = indAll([indTop(:)]);
            indLogical = logical(ones(3,N));
            indLogical(:,indRemove) = logical(0);
            indLogical = indLogical(:);
            
            
            %% find elements close to constraints and make them softer
            soft_points = find(abs(nodeM(:,3)-max(nodeM(:,3))) < 0.2);
            
            soft_points_height_threshold = min(soft_points);
            
            Y_list = Y * ones(size(elem,1),1);
            P_list = P * ones(size(elem,1),1);
            rho_list = rho * ones(size(elem,1),1);
            
            for i = 1:size(elem,1)
                if all(elem(i,:) > soft_points_height_threshold)
                    Y_list(i) = Y_list(i) / 10;
                end
            end
        case 'hete2'
            bottom_points = find(abs(nodeM(:,3)-min(nodeM(:,3))) < 0.2);
            top_points = find(abs(nodeM(:,3)-max(nodeM(:,3))) < 0.2);
            
            %% get logical indices
            indBottom = bottom_points;
            indTop = top_points;
            constraint_points = top_points;
            nFixed = length(indTop);
            
            N = size(nodeM,1); % number of nodes
            indAll = 1:N;
            indRemove = indAll([indTop(:)]);
            indLogical = logical(ones(3,N));
            indLogical(:,indRemove) = logical(0);
            indLogical = indLogical(:);
            
            
            %% find elements close to constraints and make them softer
            soft_points = find(abs(nodeM(:,3)-max(nodeM(:,3))) < 0.35);
            
            soft_points_height_threshold = min(soft_points);
            
            Y_list = Y * ones(size(elem,1),1);
            P_list = P * ones(size(elem,1),1);
            rho_list = rho * ones(size(elem,1),1);
            
            for i = 1:size(elem,1)
                if all(elem(i,:) > soft_points_height_threshold)
                    Y_list(i) = Y_list(i) / 10;
                end
            end
        case 'hete3'
            bottom_points = find(abs(nodeM(:,3)-min(nodeM(:,3))) < 0.3);
            top_points = find(abs(nodeM(:,3)-max(nodeM(:,3))) < 0.3);
            
            %% get logical indices
            indBottom = bottom_points;
            indTop = top_points;
            constraint_points = top_points;
            nFixed = length(indTop);
            
            N = size(nodeM,1); % number of nodes
            indAll = 1:N;
            indRemove = indAll([indTop(:)]);
            indLogical = logical(ones(3,N));
            indLogical(:,indRemove) = logical(0);
            indLogical = indLogical(:);
            
            
            %% find elements close to constraints and make them softer
            soft_points = find(abs(nodeM(:,3)-max(nodeM(:,3))) < 0.45);
            
            soft_points_height_threshold = min(soft_points);
            
            Y_list = Y * ones(size(elem,1),1);
            P_list = P * ones(size(elem,1),1);
            rho_list = rho * ones(size(elem,1),1);
            
            for i = 1:size(elem,1)
                if all(elem(i,:) > soft_points_height_threshold)
                    Y_list(i) = Y_list(i) / 10;
                end
            end
        case 'hete4'
            bottom_points = find(abs(nodeM(:,3)-min(nodeM(:,3))) < 0.3);
            top_points = find(abs(nodeM(:,3)-max(nodeM(:,3))) < 0.001);
            
            %% get logical indices
            indBottom = bottom_points;
            indTop = top_points;
            constraint_points = top_points;
            nFixed = length(indTop);
            
            N = size(nodeM,1); % number of nodes
            indAll = 1:N;
            indRemove = indAll([indTop(:)]);
            indLogical = logical(ones(3,N));
            indLogical(:,indRemove) = logical(0);
            indLogical = indLogical(:);
            
            
            %% find elements close to constraints and make them softer
            soft_points = find(abs(nodeM(:,3)-max(nodeM(:,3))) < 0.1);
            
            soft_points_height_threshold = min(soft_points);
            
            Y_list = Y * ones(size(elem,1),1);
            P_list = P * ones(size(elem,1),1);
            rho_list = rho * ones(size(elem,1),1);
            
            for i = 1:size(elem,1)
                if all(elem(i,:) > soft_points_height_threshold)
                    Y_list(i) = Y_list(i) / 10;
                end
            end
        case 'hete5'
            left_points = find(abs(nodeM(:,1)-min(nodeM(:,1))) < 0.003);
            right_points = find(abs(nodeM(:,1)-max(nodeM(:,1))) < 0.003);
            
            %% get logical indices
            %             indBottom = bottom_points;
            %             indTop = top_points;
            constraint_points = right_points;
            nFixed = length(indTop);
            
            N = size(nodeM,1); % number of nodes
            indAll = 1:N;
            indRemove = indAll([right_points(:)]);
            indLogical = logical(ones(3,N));
            indLogical(:,indRemove) = logical(0);
            indLogical = indLogical(:);
            
            
            %% find elements close to constraints and make them softer
            soft_points = find(abs(nodeM(:,1)-max(nodeM(:,1))) < 0.1);
            
            soft_points_height_threshold = min(soft_points);
            
            Y_list = Y * ones(size(elem,1),1);
            P_list = P * ones(size(elem,1),1);
            rho_list = rho * ones(size(elem,1),1);
            
            for i = 1:size(elem,1)
                if all(elem(i,:) > soft_points_height_threshold)
                    Y_list(i) = Y_list(i) / 10;
                end
            end
        case 'unif1'
            % constraint more points on the left
            top_points = find(abs(nodeM(:,3)-max(nodeM(:,3))) < 0.12);
            
            %% get logical indices
            %             indBottom = bottom_points;
            indTop = top_points;
            constraint_points = top_points;
            nFixed = length(indTop);
            
            N = size(nodeM,1); % number of nodes
            indAll = 1:N;
            indRemove = indAll([indTop(:)]);
            indLogical = logical(ones(3,N));
            indLogical(:,indRemove) = logical(0);
            indLogical = indLogical(:);
            
        case 'unif2'
            % constraint more points on the left
            top_points = find(abs(nodeM(:,3)-max(nodeM(:,3))) < 0.3);
            
            %% get logical indices
            %             indBottom = bottom_points;
            indTop = top_points;
            constraint_points = top_points;
            nFixed = length(indTop);
            
            N = size(nodeM,1); % number of nodes
            indAll = 1:N;
            indRemove = indAll([indTop(:)]);
            indLogical = logical(ones(3,N));
            indLogical(:,indRemove) = logical(0);
            indLogical = indLogical(:);
        case 'unif3'
            % constraint more points on the left
            top_points = find(abs(nodeM(:,3)-max(nodeM(:,3))) < 0.08);
            
            %% get logical indices
            %             indBottom = bottom_points;
            indTop = top_points;
            constraint_points = top_points;
            nFixed = length(indTop);
            
            N = size(nodeM,1); % number of nodes
            indAll = 1:N;
            indRemove = indAll([indTop(:)]);
            indLogical = logical(ones(3,N));
            indLogical(:,indRemove) = logical(0);
            indLogical = indLogical(:);
        case 'unif4'
            % constraint more points on the left
            top_points = find(abs(nodeM(:,3)-max(nodeM(:,3))) < 0.02);
            
            %% get logical indices
            %             indBottom = bottom_points;
            indTop = top_points;
            constraint_points = top_points;
            nFixed = length(indTop);
            
            N = size(nodeM,1); % number of nodes
            indAll = 1:N;
            indRemove = indAll([indTop(:)]);
            indLogical = logical(ones(3,N));
            indLogical(:,indRemove) = logical(0);
            indLogical = indLogical(:);
        case 'unif5'
            % constraint more points on the left
            top_points = find(abs(nodeM(:,3)-max(nodeM(:,3))) < 0.01);
            
            %% get logical indices
            %             indBottom = bottom_points;
            indTop = top_points;
            constraint_points = top_points;
            nFixed = length(indTop);
            
            N = size(nodeM,1); % number of nodes
            indAll = 1:N;
            indRemove = indAll([indTop(:)]);
            indLogical = logical(ones(3,N));
            indLogical(:,indRemove) = logical(0);
            indLogical = indLogical(:);
        case 'unif6'
            % constraint more points on the left
            top_points = find(abs(nodeM(:,3)-max(nodeM(:,3))) < 0.3);
            
            %% get logical indices
            %             indBottom = bottom_points;
            indTop = top_points;
            constraint_points = top_points;
            nFixed = length(indTop);
            
            N = size(nodeM,1); % number of nodes
            indAll = 1:N;
            indRemove = indAll([indTop(:)]);
            indLogical = logical(ones(3,N));
            indLogical(:,indRemove) = logical(0);
            indLogical = indLogical(:);
    end
    
elseif any(strfind(meshfile,'flatTorusmesh'))
    switch scene_name
        case 'hete1'
            left_points = find(abs(nodeM(:,1)-min(nodeM(:,1))) < 0.003);
            right_points = find(abs(nodeM(:,1)-max(nodeM(:,1))) < 0.003);
            
            %% get logical indices
            %             indBottom = bottom_points;
            %             indTop = top_points;
            constraint_points = right_points;
            
            
            N = size(nodeM,1); % number of nodes
            indAll = 1:N;
            indRemove = indAll([right_points(:)]);
            indLogical = logical(ones(3,N));
            indLogical(:,indRemove) = logical(0);
            indLogical = indLogical(:);
            
            
            %% find elements close to constraints and make them softer
            soft_points = find(abs(nodeM(:,1)-max(nodeM(:,1))) < 0.1);
            
            soft_points_height_threshold = min(soft_points);
            
            Y_list = Y * ones(size(elem,1),1);
            P_list = P * ones(size(elem,1),1);
            rho_list = rho * ones(size(elem,1),1);
            
            for i = 1:size(elem,1)
                if all(elem(i,:) > soft_points_height_threshold)
                    Y_list(i) = Y_list(i) / 10;
                end
            end
            case 'hete2'
            left_points = find(abs(nodeM(:,1)-min(nodeM(:,1))) < 0.1);
            right_points = find(abs(nodeM(:,1)-max(nodeM(:,1))) < 0.1);
            
            %% get logical indices
            %             indBottom = bottom_points;
            %             indTop = top_points;
            constraint_points = right_points;
            
            
            N = size(nodeM,1); % number of nodes
            indAll = 1:N;
            indRemove = indAll([right_points(:)]);
            indLogical = logical(ones(3,N));
            indLogical(:,indRemove) = logical(0);
            indLogical = indLogical(:);
            
            
            %% find elements close to constraints and make them softer
            soft_points = find(abs(nodeM(:,1)-max(nodeM(:,1))) < 0.2);
            
            soft_points_height_threshold = min(soft_points);
            
            Y_list = Y * ones(size(elem,1),1);
            P_list = P * ones(size(elem,1),1);
            rho_list = rho * ones(size(elem,1),1);
            
            for i = 1:size(elem,1)
                if all(elem(i,:) > soft_points_height_threshold)
                    Y_list(i) = Y_list(i) / 10;
                end
            end
        case 'unif1'
            % constraint more points on the left
            top_points = find(abs(nodeM(:,1)-max(nodeM(:,1))) < 0.12);
            
            %% get logical indices
            %             indBottom = bottom_points;
            indTop = top_points;
            constraint_points = top_points;
            nFixed = length(indTop);
            
            N = size(nodeM,1); % number of nodes
            indAll = 1:N;
            indRemove = indAll([indTop(:)]);
            indLogical = logical(ones(3,N));
            indLogical(:,indRemove) = logical(0);
            indLogical = indLogical(:);
    end
end

end
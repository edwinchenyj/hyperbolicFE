function test_modes_free()
clear all
close all


nodeM = [                % list of xy "node" coordinates
    1, 0
    1, 1
    0, 1
    0, 0
    ] ;

edge = [                % list of "edges" between nodes
    1, 2                % outer square
    2, 3
    3, 4
    4, 1
    ] ;

hfun_list = [1 0.5 0.2 0.1 0.05 0.02]'; % list of target edge length, the parameter to change mesh resolution
sm_eig_list = zeros(length(hfun_list),4); % allocate the space for the smallest eigenvalues
total_mass_list = zeros(length(hfun_list),1);
for i = 1:length(hfun_list)

    hfun = hfun_list(i) ; % target edge length, the parameter to change mesh resolution
[nodeM,etri, ...
    elem,tnum] = refine2(nodeM,edge,[],[],hfun) ;

% swap the elem order to make area positive (the convention used happened
% to be different
elem(:,[1 3]) = elem(:,[3 1]);

N = size(nodeM,1);

% pin all the points on line x = -1
% construct the logical index for the free points
Xind_fix = nodeM(:,1) == -1; % logical indices for x
ind_fix = reshape(transpose(repmat(Xind_fix,1,2)),[],1); % logical index for total position vector


% construct triangular mesh object
obj = staticTriMesh(nodeM, elem);

Y = 1; % Young's modululs
P = 0.4; % Poisson ratio
rho = 1; % density

obj.SetMaterial( Y, P, rho, 1:size(elem,1), 1); % set the tets to be neohookean

M = obj.M;
% print out the total mass for sanity check
% should be the same for all mesh resolution
% disp('total mass')
% disp(sum(spdiags(M)))
total_mass_list(i) = sum(spdiags(M));

M = M(~ind_fix,~ind_fix); % extract the non-fixed part

Dx = 0*rand(2*N,1); % displacement field. set zero for rest state
obj.SetCurrentState(Dx);

K = obj.StiffnessMatrix;
K = K(~ind_fix,~ind_fix); % extract the non-fixed part

% disp('free mass')
% disp(sum(spdiags(M)))

% extract the 4 smallest eigenvalues and eigenvectors
[v,d] = eigs(K,4,'sm');
% disp('smallest eigs')
sm_eig = diag(d);
% disp(sort(sm_eig))
sm_eig_list(i,:) = sort(sm_eig)';

% extract the smallest non-zero eigenvalue and eigenvector
sm_eig(~(abs(sm_eig)>1e-4)) = 0;
ind = sm_eig == min(sm_eig(sm_eig>0));
sm_eigv = v(:,ind);

% use the eigenvector for the heightmap and coloring
% the fixed nodes are 0
figure('Position',[500+10*i,100-10*i,1000,800])
subplot(2,2,1)
trimesh(elem,nodeM(:,1),nodeM(:,2))
title('Mesh top view')
axis equal

subplot(2,2,3)
Cx = zeros(length(nodeM(:,1)),1);
Cx(~Xind_fix) = sm_eigv(1:2:end); % x component of the eigenvector
trimesh(elem,nodeM(:,1),nodeM(:,2),Cx,Cx,'FaceColor','interp');
title('lowest eigen mode X component')
dim = [0.5 0.1 0.3 0.3];
str = {'hfun =',num2str(hfun), 'eig = ', num2str(sm_eig(ind))};
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',14,'Color','k');
axis equal

subplot(2,2,4)
Cy = zeros(length(nodeM(:,1)),1);
Cy(~Xind_fix) = sm_eigv(2:2:end); % y component of the eigenvector
trimesh(elem,nodeM(:,1),nodeM(:,2),Cy,Cy,'FaceColor','interp');
title('lowest eigen mode Y component')
dim = [0.8 0.1 0.3 0.3];
str = {'hfun =',num2str(hfun), 'eig = ', num2str(sm_eig(ind))};
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',14,'Color','k');
% trimesh(elem,nodeM(:,1),nodeM(:,2),zeros(size(nodeM,1),1),sm_eigv(1:2:end));
axis equal

subplot(2,2,2)
trimesh(elem,nodeM(:,1) + Cx ,nodeM(:,2) + Cy)
title('Mesh top view')
axis equal
drawnow

obj.mesh_quality();
end
eval(['table(hfun_list, total_mass_list,', sprintf('sm_eig_list(:,%d),',1:length(sm_eig)), '''VariableNames'', {''EdgeLength'' ''TotalMass'' ', sprintf('''Eigs%d'' ',1:length(sm_eig)), '})'])

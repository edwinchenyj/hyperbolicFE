% using distmesh to generate tet meshes and 2D surface meshes

%% ball

type = 'ball';
radius = 0.2;
n = 5; % the parameter that controls how fine the mesh is

fs = filesep;
meshname = sprintf('simulation data%cmeshes%c%s r%.e n%d',fs,fs,type, radius,n);

if exist([meshname '.mat'], 'file') ~= 2
rng('default'); % Always the same results
rng(1);
fd=@(p) dsphere(p,0,0,0,1);
[p,t]=distmeshnd(fd,@huniform,1/n,[-1,-1,-1;1,1,1],[]);
[p,t] = fixmesh(p,t); % fix the orientation
nodeM = radius * p;
elem = t;
else
    load([meshname '.mat'], 'nodeM', 'elem');
end
colormap jet
s_alpha = sprintf('%f',1/n);
% figure
tetramesh(elem,nodeM,ones(size(elem,1)),'FaceAlpha',s_alpha,'EdgeAlpha',s_alpha);
axis equal
grid on
% figure
% simpplot(nodeM,elem)
save([meshname '.mat'],'nodeM','elem')

%% torus mesh
% 
type = 'torus';
radius = 5;
n = 6; % the parameter that controls how detail the mesh is
h = 1/n;
fs = filesep;
meshname = sprintf('%s r%d n%d',type, radius,n);

if exist([meshname '.mat'], 'file') ~= 2
rng('default'); % Always the same results
rng(1);
fd=@(p) (sum(p.^2,2)+.8^2-.2^2).^2-4*.8^2*(p(:,1).^2+p(:,2).^2);
    [p,t]=distmeshsurface(fd,@huniform,h,[-1.1,-1.1,-.25;1.1,1.1,.25]);
    [p,t] = fixmesh(p,t); % fix the orientation
nodeM = radius * p;
elem = t;
else
    load([meshname '.mat'], 'nodeM', 'elem');
end
simpplot(nodeM,elem);
axis equal
grid on
save([meshname '.mat'],'nodeM','elem')

%% Tet bar from distmesh
type = 'bar';
h = 0.2;
leng = 1;
fs = filesep;
meshname = sprintf('simulation data%cmeshes%c%s r%d n%d',fs,fs,type, h, leng);

box = [
    0.0, 0.0, 0.0;
    3.0, 1.0, 1.0 ];


iteration_max = 200;

pfix = [
    0.0, 0.0, 0.0;
    0.0, 0.0, 1.0;
    0.0, 1.0, 0.0;
    0.0, 1.0, 1.0;
    3.0, 0.0, 0.0;
    3.0, 0.0, 1.0;
    3.0, 1.0, 0.0;
    3.0, 1.0, 1.0 ];

[ p, t ] = distmesh_3d ( @fd01, @fh01, h, box, iteration_max, pfix );

[p,t] = fixmesh(p,t); % fix the orientation
nodeM = p /3 /5 * leng;
elem = t;

simpplot(nodeM,elem,[],[],[],[],[],1);
axis equal
grid on
save([meshname '.mat'],'nodeM','elem')
%% Tet large bar from distmesh
type = 'large_bar';
h = 0.2;
leng = 1;
fs = filesep;
meshname = sprintf('simulation data%cmeshes%c%s r%d n%d',fs,fs,type, h, leng);

box = [
    0.0, 0.0, 0.0;
    3.0, 1.0, 1.0 ];


iteration_max = 200;

pfix = [
    0.0, 0.0, 0.0;
    0.0, 0.0, 1.0;
    0.0, 1.0, 0.0;
    0.0, 1.0, 1.0;
    3.0, 0.0, 0.0;
    3.0, 0.0, 1.0;
    3.0, 1.0, 0.0;
    3.0, 1.0, 1.0 ];

[ p, t ] = distmesh_3d ( @fd01, @fh01, h, box, iteration_max, pfix );

[p,t] = fixmesh(p,t); % fix the orientation
nodeM = p * 10 * leng;
elem = t;

simpplot(nodeM,elem);
axis equal
grid on
save([meshname '.mat'],'nodeM','elem')
%% Tet long bar from distmesh
type = 'long_bar';
h = 0.2;
leng = 1;
fs = filesep;
meshname = sprintf('simulation data%cmeshes%c%s r%d n%d',fs,fs,type, h, leng);

box = [
    0.0, 0.0, 0.0;
    5.0, 1.0, 1.0 ];


iteration_max = 200;

pfix = [
    0.0, 0.0, 0.0;
    0.0, 0.0, 1.0;
    0.0, 1.0, 0.0;
    0.0, 1.0, 1.0;
    5.0, 0.0, 0.0;
    5.0, 0.0, 1.0;
    5.0, 1.0, 0.0;
    5.0, 1.0, 1.0 ];

[ p, t ] = distmesh_3d ( @fd01_long, @fh01, h, box, iteration_max, pfix );

[p,t] = fixmesh(p,t); % fix the orientation
nodeM = p /5 /5 * leng;
elem = t;

simpplot(nodeM,elem);
axis equal
grid on
save([meshname '.mat'],'nodeM','elem')
%% Tet small bar from distmesh
type = 'small_bar';
h = 0.1;
leng = 1;
fs = filesep;
meshname = sprintf('simulation data%cmeshes%c%s r%d n%d',fs,fs,type, h, leng);

box = [
    0.0, 0.0, 0.0;
    3.0, 1.0, 1.0 ];


iteration_max = 200;

pfix = [
    0.0, 0.0, 0.0;
    0.0, 0.0, 1.0;
    0.0, 1.0, 0.0;
    0.0, 1.0, 1.0;
    3.0, 0.0, 0.0;
    3.0, 0.0, 1.0;
    3.0, 1.0, 0.0;
    3.0, 1.0, 1.0 ];

[ p, t ] = distmesh_3d ( @fd01, @fh01, h, box, iteration_max, pfix );

[p,t] = fixmesh(p,t); % fix the orientation
nodeM = p /3 * leng / 5;
elem = t;

simpplot(nodeM,elem);
axis equal
grid on
save([meshname '.mat'],'nodeM','elem')

%% Tet cube from distmesh
type = 'cube';
h = 0.2;
leng = 1;
fs = filesep;
meshname = sprintf('simulation data%cmeshes%c%s r%d n%d',fs,fs,type, h, leng);
  
  box = [ 
      0.0, 0.0, 0.0;
      1.0, 1.0, 1.0 ];
  
iteration_max = 200;

  pfix = [
      0.0, 0.0, 0.0;
      0.0, 0.0, 1.0;
      0.0, 1.0, 0.0;
      0.0, 1.0, 1.0;
      1.0, 0.0, 0.0;
      1.0, 0.0, 1.0;
      1.0, 1.0, 0.0;
      1.0, 1.0, 1.0 ];

  [ p, t ] = distmesh_3d ( @fd03, @fh03, h, box, iteration_max, pfix );

[p,t] = fixmesh(p,t); % fix the orientation
nodeM = p * leng;
elem = t;

simpplot(nodeM,elem);
axis equal
grid on
save([meshname '.mat'],'nodeM','elem')

%% Tet cylinder from distmesh
type = 'cylinder';
h = 0.2;
leng = 1;
fs = filesep;
meshname = sprintf('simulation data%cmeshes%c%s r%d n%d',fs,fs,type, h, leng);
  
 
  box = [ 
      0.0, 0.0, 0.0;
      1.0, 1.0, 4.0 ];
  
iteration_max = 800;

  pfix = [];

  [ p, t ] = distmesh_3d ( @fd02, @fh02, h, box, iteration_max, pfix );
  
[p,t] = fixmesh(p,t); % fix the orientation
nodeM = p /4 * leng;
elem = t;

simpplot(nodeM,elem);
axis equal
grid on
save([meshname '.mat'],'nodeM','elem')

%% Tetrahedron bar
type = 'bar';
leng = 1;
width = 0.1;
for elem_size = [0.1,0.05,0.02,0.01];
    
    fs = filesep;
    meshname = sprintf('simulation data%cmeshes%c%s l%d es%d',fs,fs,type, leng,100*elem_size);
    
    if exist([meshname '.mat'], 'file') ~= 2
        
        
        xd = 0:elem_size:width;
        yd = 0:elem_size:leng;
        zd = 0:elem_size:width;
        [x,y,z] = meshgrid(xd,yd,zd); % a cube
        x = [x(:)];
        y = [y(:)];
        z = [z(:)];
        nodeM = [x,y,z];
        elem = delaunayTriangulation(x,y,z);
        
    else
        load([meshname '.mat'], 'nodeM', 'elem');
    end
    simpplot(nodeM,elem);
    axis equal
    grid on
    save(meshname,'nodeM','elem')
    1
end

%% Tetrahedron bar
type = 'bar';
leng = 2;
width = 0.1;
for elem_size = [0.1,0.05,0.02,0.01];
    
    fs = filesep;
    meshname = sprintf('simulation data%cmeshes%c%s l%d es%d',fs,fs,type, leng,100*elem_size);
    
    if exist([meshname '.mat'], 'file') ~= 2
        
        
        xd = 0:elem_size:width;
        yd = 0:elem_size:leng;
        zd = 0:elem_size:width;
        [x,y,z] = meshgrid(xd,yd,zd); % a cube
        x = [x(:)];
        y = [y(:)];
        z = [z(:)];
        nodeM = [x,y,z];
        elem = delaunayTriangulation(x,y,z);
        
    else
        load([meshname '.mat'], 'nodeM', 'elem');
    end
    simpplot(nodeM,elem);
    axis equal
    grid on
    save(meshname,'nodeM','elem')
    1
end

%%
fd=@(p) (sum(p.^2,2)+.8^2-.2^2).^2-4*.8^2*(p(:,1).^2+p(:,2).^2);
    [p,t]=distmeshsurface(fd,@huniform,0.1,[-1.1,-1.1,-.25;1.1,1.1,.25]);
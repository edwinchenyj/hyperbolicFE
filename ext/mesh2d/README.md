# `MESH2D: A MATLAB-based mesh generator`

`MESH2D` is a simple `MATLAB`-based Delaunay mesh-generator for two-dimensional geometries. It is designed to generate high-quality constrained Delaunay triangulations for general polygonal regions in the plane. `MESH2D` provides simple and yet effective implementations of "Delaunay-refinement" and "Frontal-Delaunay" triangulation techniques, in additional to "hill-climbing" type mesh-optimisation. Support for user-defined "mesh-spacing" functions and "multi-part" geometry definitions are provided, allowing varying levels of mesh-resolution to be specified within complex domains.

Algorithms implemented in `MESH2D` are "provably-good" - ensuring convergence, geometrical and topological correctness, and providing guarantees on worst-case element quality in all cases.

<p align="center">
  <img src = "../master/poly-data/lake-1-small.png"> &nbsp &nbsp &nbsp &nbsp
  <img src = "../master/poly-data/lake-2-small.png">
</p>

The `MESH2D` package has been around now for a while(!), having been first released in 2006 as an extension of the <a href="http://persson.berkeley.edu/distmesh/">`DISTMESH`</a> algorithm of Persson and Strang. `MESH2D-3.0.0` is a major re-release, representing a complete re-working of the code to (a) bring it up-to-date with more recent developments in meshing research, and (b) utilise the robust constrained Delaunay triangulation support now provided by `MATLAB`. `MESH2D-3.0.0` is essentially now a simplified version of my <a href="https://github.com/dengwirda/jigsaw-matlab/">`JIGSAW`</a> mesh-generation algorithm, providing an implementation of "provably-good" Delaunay-refinement and Frontal-Delaunay triangulation techniques.

# `Starting Out`

After downloading and unzipping the current <a href="https://github.com/dengwirda/mesh2d/archive/master.zip">repository</a>, navigate to the installation directory within <a href="http://www.mathworks.com">`MATLAB`</a> / <a href="https://www.gnu.org/software/octave">`OCTAVE`</a> and run the set of examples contained in `tridemo.m`:
````
tridemo( 0); % a very simple example to get everything started.
tridemo( 1); % investigate the impact of the "radius-edge" threshold.
tridemo( 2); % Frontal-Delaunay vs. Delaunay-refinement refinement.
tridemo( 3); % explore impact of user-defined mesh-size constraints.
tridemo( 4); % explore impact of "hill-climbing" mesh optimisations.
tridemo( 5); % assemble triangulations for "multi-part" geometries.
tridemo( 6); % assemble triangulations with "internal" constraints.
tridemo( 7); % investigate the use of "quadtree"-type refinement.
tridemo( 8); % explore use of custom, user-defined mesh-size functions.
tridemo( 9); % larger-scale problem, mesh refinement + optimisation. 
tridemo(10); % medium-scale problem, mesh refinement + optimisation. 
````
# `Attribution!`

If you make use of `MESH2D` please include a reference to the following! `MESH2D` is designed to provide a simple and easy-to-understand implementation of Delaunay-based mesh-generation techniques. For a much more advanced, and fully three-dimensional mesh-generation library, see the <a href="https://github.com/dengwirda/jigsaw-matlab/">`JIGSAW`</a> package. `MESH2D` makes use of the <a href="https://github.com/dengwirda/aabb-tree">`AABBTREE`</a> and <a href="https://github.com/dengwirda/find-tria">`FINDTRIA`</a> packages to compute efficient spatial queries and intersection tests. 

`[1]` - Darren Engwirda, <a href="http://hdl.handle.net/2123/13148">Locally-optimal Delaunay-refinement and optimisation-based mesh generation</a>, Ph.D. Thesis, School of Mathematics and Statistics, The University of Sydney, September 2014.

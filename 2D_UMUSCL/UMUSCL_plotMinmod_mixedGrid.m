%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%        Unstructured MUSCL solver for 2-d scalar advection equation
%                      by Manuel Diaz, NTU, 29.04.2015
%
%                          u_t + u_x + u_y = 0,
%
%   MUSCL based numerical schemes extend the idea of using a linear
%   piecewise approximation to each cell by using slope limited left and
%   right extrapolated states. This results in the following high
%   resolution, TVD discretisation scheme.   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Refs:
%   [1] Wikipedia, MUSCL scheme, available online at:
%   http://en.wikipedia.org/wiki/MUSCL_scheme
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;

% Build grid of mixed cells
[vx,vy,EtoV,nE,nN,BC] = SquareMesh_v2('MIXED');

% Build Unstructured Mesh
[node,elem,edge] = BuildUnstructuredMesh2d(vx,vy,EtoV,nE,nN,BC);

% Build Cell Center's grid
Xc = zeros(nE,1); Yc = zeros(nE,1);
for e = 1:nE
    Xc(e) = elem(e).x; Yc(e) = elem(e).y;
end

% Build Cell Average's IC 
IC = 02;
switch IC
    case 01 % square jump
        f0 = @(x,y) zeros(size(x)) + 1*(x>=0.3 & x<0.7).*(y>=0.3 & y<0.7);
    case 02 % gaussian
        f0 = @(x,y) exp(-30*((x-1/2).^2+(y-1/2).^2));
end
u0 = f0(Xc,Yc);

% Visualize Initial Condition
figure(2); hold on;
for e = 1:nE
    patch(vx(EtoV{e})',vy(EtoV{e})',f0(vx(EtoV{e})',vy(EtoV{e})'),u0(e));
end
grid on; hold off

%Visualize Initial Condition as cell averages
figure(3); hold on;
for e = 1:nE
    %u0(e) = f0(Xc(e),Yc(e));
    u0(e) = cellAvg2(f0,[vx(EtoV{e}),vy(EtoV{e})]);
    patch(vx(EtoV{e})',vy(EtoV{e})',u0(e)*ones(1,numel(EtoV{e})),u0(e));
end
grid on; hold off

% Build Structure Arrays
cell(nE).u=0; 
for e = 1:nE % for every cell
    cell(e).x = Xc(e);
    cell(e).y = Yc(e);
    cell(e).u = u0(e);
    cell(e).vx= [node(elem(e).v).x];
    cell(e).vy= [node(elem(e).v).y];
    cell(e).ui= zeros(elem(e).nV,1);
end

% Build Left and Right slopes for every element
bCells = zeros(nE,1);   % list of boundary cells
for ei = 1:nE
    ni = elem(ei).nV;                   % element total edges/neighbors
    ki = elem(ei).nghbr; ki=ki(ki~=0);  % non-zero neighbours elements
    if  numel(ki) == ni                 % only non-boundary elements
        % cells center coordiantes and cell average
        xi = elem(ei).x;
        yi = elem(ei).y;
        ui = cell(ei).u;
        % neighbours cell center coordinates and cell averages
        xk = [elem(ki).x];
        yk = [elem(ki).y];
        uk = [cell(ki).u];
        % build inv(AtA) matrix
        invAtA = LSQinvMat2d(xi,yi,ei,xk,yk);
        gradu = LSQgradients2d(xi,yi,ui,xk,yk,uk,invAtA);
        cell(ei).ux= gradu(1);
        cell(ei).uy= gradu(2);
    else
        bCells(ei) = 1; % Mark as boundary cell
    end
end
bCells = find(bCells);

% Build Left and Right extrapolations for every cell
for e = setdiff(1:nE,bCells)
    for v = 1:elem(e).nV
       cell(e).ui(v) = cell(e).u + ...
           cell(e).ux*(cell(e).vx(v)-cell(e).x)+ ...
           cell(e).uy*(cell(e).vy(v)-cell(e).y);
    end
end

% Plot Piezewise Linear Approximation with limiter
% Visualize Initial Condition
figure(4); hold on;
for e = setdiff(1:nE,bCells)
    patch(cell(e).vx,cell(e).vy,cell(e).ui,u0(e));
end
grid on; hold off

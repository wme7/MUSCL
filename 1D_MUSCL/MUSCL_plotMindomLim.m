%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%               basic MUSCL solver for scalar advection equation
%                      by Manuel Diaz, NTU, 29.04.2015
%
%                             u_t + u_x = 0,
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

% Set number of cells
nx = 21; 

% Build Cell Center's grid
lx=1; dx=lx/nx; xc=(dx/2:dx:1);

% Build Cell Average's IC 
IC = 01;
switch IC
    case 01 % square jump
        u0 = zeros(size(xc)) + 1*(xc>=0.3 & xc<0.7); range=[0,1,-0.25,1.25];
    case 02 % gaussian
        u0 = exp(-30*(xc-0.5).^2); range=[0,1,-0.25,1.25];
end

% Visualize Initial Condition
plot(xc,u0,'-k','LineWidth',1.5); hold on;

% Build Structure Arrays
cell(nx).BC=1; cell(1).BC=1;
for j = 1:nx % for every cell
    cell(j).x = xc(j);
    cell(j).u = u0(j);
end

for j = 2:nx-1 % internal cells
    % Left and Right slopes 
    cell(j).dul = (cell( j ).u - cell(j-1).u)/dx; % simple diff between j and j-1
    cell(j).dur = (cell(j+1).u - cell( j ).u)/dx; % simple diff between j+1 and j
    cell(j).du = minmod([cell(j).dul,cell(j).dur]);
    %cell(j).du = mean([cell(j).dul,cell(j).dur]);
    cell(j).dx = dx;
end

%     j+1/2         Cell's grid:
%   | wL|   |
%   |  /|wR |           1   2   3   4        N-2 N-1  N
%   | / |\  |   {x=0} |-o-|-o-|-o-|-o-| ... |-o-|-o-|-o-| {x=L}
%   |/  | \ |         1   2   3   4   5        N-1  N  N+1
%   |   |  \|
%   |   |   |       NC: Here cells 1 and N are ghost cells!
%     j  j+1
%
% Build Left and Right extrapolations for every cell
for j = 2:nx-1 % internal cells
    cell(j).uR = cell(j).u - cell(j).du*dx/2; % u_{j-1/2}^{+} from j 
    cell(j).uL = cell(j).u + cell(j).du*dx/2; % u_{j+1/2}^{-} from j
end
disp(struct2table(cell));

% Plot Piezewise Linear Approximation with limiter
xi= [-1/2;0;1/2];
x = repmat([cell(2:nx-1).x],3,1) + xi*[cell.dx];	
u = [cell(2:nx-1).uR;cell(2:nx-1).u;cell(2:nx-1).uL];
plot(x,u,'-o','LineWidth',1.5); grid on; axis(range); hold off
title('Testing the minmod function in 1d');
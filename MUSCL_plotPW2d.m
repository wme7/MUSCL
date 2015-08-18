%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%             basic MUSCL solver for 2-d scalar advection equation
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

% Set number of cells
nx = 20; ny = 20;

% Build Cell Center's grid
lx=1; dx=lx/nx; xc=(dx/2:dx:1);
ly=1; dy=ly/ny; yc=(dy/2:dy:1);
[Xc,Yc] = meshgrid(xc,yc);

% Build Cell Average's IC 
IC = 02;
switch IC
    case 01 % square jump
        u0 = zeros(size(Xc)) + 1*(Xc>=0.3 & Xc<0.7).*(Yc>=0.3 & Yc<0.7);
    case 02 % gaussian
        u0 = exp(-30*((Xc-1/2).^2+(Yc-1/2).^2));
end

% Visualize Initial Condition
subplot(121); surf(Xc,Yc,u0); view(-20,20);

% Build Structure Arrays
cell(nx,ny).BC=1; cell(1,1).BC=1;
for i = 1:ny; % for every cell
    for j = 1:nx; % for every cell
        cell(i,j).x = Xc(i,j);
        cell(i,j).y = Yc(i,j);
        cell(i,j).u = u0(i,j);
    end
end

for i = 2:ny-1
    for j = 2:nx-1 % internal cells
        % Left and Right slopes 
        cell(i,j).duw = (cell(i, j ).u - cell(i,j-1).u)/dx; % du between j and j-1
        cell(i,j).due = (cell(i,j+1).u - cell(i, j ).u)/dx; % du between j+1 and j
        %cell(i,j).dux = minmod([cell(i,j).duw,cell(i,j).due]);
        cell(i,j).dux = mean([cell(i,j).duw,cell(i,j).due]);
        cell(i,j).dus = (cell( i ,j).u - cell(i-1,j).u)/dy; % du between i and i-1
        cell(i,j).dun = (cell(i+1,j).u - cell( i ,j).u)/dy; % du between i+1 and i
        %cell(i,j).duy = minmod([cell(i,j).dus,cell(i,j).dun]);
        cell(i,j).duy = mean([cell(i,j).dus,cell(i,j).dun]);
        cell(i,j).dy = dy; cell(i,j).dx = dx;
    end
end

% 1-D perspective:
%
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
for i = 2:ny-1
    for j = 2:nx-1 % internal cells
     % extrapolation to the sides
     cell(i,j).uW = cell(i,j).u - cell(i,j).dux*dx/2; % u_{i,j-1/2}^{-} from (i,j) 
     cell(i,j).uE = cell(i,j).u + cell(i,j).dux*dx/2; % u_{i,j+1/2}^{-} from (i,j)
     cell(i,j).uS = cell(i,j).u - cell(i,j).duy*dy/2; % u_{i-1/2,j}^{-} from (i,j) 
     cell(i,j).uN = cell(i,j).u + cell(i,j).duy*dy/2; % u_{i+1/2,j}^{-} from (i,j)
     % extrapolation to the corners
     cell(i,j).uNE = cell(i,j).u + cell(i,j).dux*dx/2 + cell(i,j).duy*dy/2; % u_{i+1/2,j+1/2}^{-} from (i,j) 
     cell(i,j).uSE = cell(i,j).u + cell(i,j).dux*dx/2 - cell(i,j).duy*dy/2; % u_{i-1/2,j+1/2}^{-} from (i,j)
     cell(i,j).uNW = cell(i,j).u - cell(i,j).dux*dx/2 + cell(i,j).duy*dy/2; % u_{i-1/2,j-1/2}^{-} from (i,j) 
     cell(i,j).uSW = cell(i,j).u - cell(i,j).dux*dx/2 - cell(i,j).duy*dy/2; % u_{i+1/2,j-1/2}^{-} from (i,j)
    end
end

% Plot Piezewise Linear Approximation with limiter
xi= [-1/2;1/2]; eta=[-1/2;1/2]; [XI,ETA] = meshgrid(xi,eta);
subplot(122); hold on; view(-20,20); grid on;
for i = 2:ny-1
    for j = 2:nx-1 % internal cells
        surf(XI*dx+cell(i,j).x,ETA*dy+cell(i,j).y,...
            [cell(i,j).uSW,cell(i,j).uSE;cell(i,j).uNW,cell(i,j).uNE])
    end
end
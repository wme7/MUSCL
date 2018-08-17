function ubar = cellAvg(f,x,y)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Compute the cell average of a triangles and quadrilateral elements
%                    Coded by Manuel Diaz, NTU 2015.05.29
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following implementation is base on the following discussions
%
% Refs: 
% [1] http://www.mathworks.com/matlabcentral/newsreader/view_thread/276943
% [2] http://stackoverflow.com/questions/14696525/double-integration
%     -over-a-polygon-in-matlab
% [3] http://connor-johnson.com/2014/03/09/integrals-over-arbitrary
%     -triangular-regions-for-fem/
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;

vx = [0,0,1];
vy = [0,1,0];
plot(vx,vy,'r-'); axis equal;

% Obtain number of vertices
p = numel(vx);

%x1 = vx(1); y1 = vy(1);
%x2 = vx(2); y2 = vy(2);
%x3 = vx(3); y3 = vy(3);

% produce integration
ubar = quad2d(fun,0,1,0,1);

% What was the triangular region in x,y space then becomes a unit square
% in u,v space - the variables u and v each range from 0 to 1.
 end

function f = fun(u,v,VX,VY)
% Define Integrand function

% Obtain vertices coordinates (they are like constants here)
x1 = VX(1); y1 = VY(1);
x2 = VX(2); y2 = VY(2);
x3 = VX(3); y3 = VY(3);

% transformation form unit square [0,1]x[0,1] to the unit triangle
x = (1-u)*x1 + u.*((1-v)*x2 + v*x3); 
y = (1-u)*y1 + u.*((1-v)*y2 + v*y3);

% The determinant of the jacobian transformation becomes,
dxdu = ( (1-v)*x2 + v*x3 - x1 );
dxdv = ( u*x3 - u*x2 );
dydu = ( (1-v)*y2 + v*y3 - y1 );
dydv = ( u*y3 - u*y2 );
detJ = ( dxdu.*dydv - dxdv.*dydu ); 

% Fucntion to be evaluated
f = x.*y./detJ;
end
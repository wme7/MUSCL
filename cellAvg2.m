% Refs:
% [1] https://www.math.ohiou.edu/courses/math3600/lecture25.pdf

function I = cellAvg2(f,V) 
% Integrates a function based on a triangulation , using three corners
% Inputs : f -- the function to integrate , as an inline
% V -- the vertices .
% Each row has the x and y coordinates of a vertex
% T -- the triangulation .
% Each row gives the indices of three corners
% Output : the approximate integral

x = V (:,1); % extract x and y coordinates of all nodes
y = V (:,2);
nv =size(V,1);
switch nv
    case 3
        %A = polyarea(x,y); % find coordinates and area
        z1 = f (x(1) , y(1) ); % find values and average
        z2 = f (x(2) , y(2) );
        z3 = f (x(3) , y(3) );
        I = ( z1 + z2 + z3 )/3;
    case 4
        %A = polyarea(x,y); % find coordinates and area
        z1 = f (x(1) , y(1) ); % find values and average
        z2 = f (x(2) , y(2) );
        z3 = f (x(3) , y(3) );
        z4 = f (x(4) , y(4) );
        I = ( z1 + z2 + z3 + z4 )/4;
end
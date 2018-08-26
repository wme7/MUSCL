function [res,wsn]=UMUSCL_AdvRHS(u,F,dF,G,dG,node,edge,bound,limiter,fluxMth)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computes the residual for a node-centered finite-volume method
% -------------------------------------------------------------------------
%  Input: the current solution
% Output: node(:)%res = the residual computed by the current solution.
% -------------------------------------------------------------------------
% Note: dU/dt + dF/dx + dG/dy = 0. Residuals are first computed as
%       the integral of (dF/dx + dG/dy), and at the end negative sign is
%       added so that we have dU/dt = Res at every node.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 0. Load u into the nodes data
for i=1:nN; node(i).u=u(i); end

% 1. Initialize residual and fluxes arrays in nodes
for i=1:nN; node(i).res=0; node(i).wsn=0; end 

% 2. Gradient Reconstruction
for i=1:nN; node(i).gradu=LSQgradients2d(node(i).x,node(i).y,node(i).w,...
        [node(node(i).nghbr).x],[node(node(i).nghbr).y],...
        [node(node(i).nghbr).w], node(i).invAtA ); 
end

% Flux computation across internal edges (to be accumulated in res(:))
%
%   node2              1. Extrapolate the solutions to the edge-midpoint
%       o                 from the nodes, n1 and n2.
%        \   face      2. Compute the numerical flux
%         \ -------c2  3. Add it to the residual for n1, and subtract it 
%        / \              from the residual for n2.
%   face/   \ edge
%      /     o         Directed area is the sum of the left and the right 
%    c1    node1       faces. Left/right face is defined by the edge-midpoint 
%                      and the centroid of the left/right element.
%                      Directed area is positive in n1 -> n2
%
% (c1, c2: element centroids)

% Left and right nodes of the i-th edge

    node1 = edge(i).n1;  % Left node of the edge
    node2 = edge(i).n2;  % Right node of the edge
      n12 = edge(i).dav; % This is the directed area vector (unit vector)
  mag_n12 = edge(i).da;  % Magnitude of the directed area vector
      e12 = edge(i).ev;  % This is the vector along the edge (uniti vector)
  mag_e12 = edge(i).e;   % Magnitude of the edge vector (Length of the edge)

% Solution gradient projected along the edge
%
%  NOTE: The gradient is multiplied by the distance.
%        So, it is equivalent to the solution difference.

duL = (node(node1).gradu(:,1)*e12(1) + node(node1).gradu(:,2)*e12(2))*0.5*mag_e12;
duR = (node(node2).gradu(:,1)*e12(1) + node(node2).gradu(:,2)*e12(2))*0.5*mag_e12;

%  It is now limiter time
%--------------------------------------------------------------------------
% In 1D: dwp = w_{j+1}-w_j, dwm = w_j-w_{j-1} 
%               => limited_slope = limiter(dwm,dwp)
%
% We can do the same in 2D as follows.
%
% In 2D:    dwp = w_{neighbor}-w_j, dwm = 2*(grad_w_j*edge)-dwp
%               => limited_slope = limiter(dwm,dwp)
%
% NOTE: On a regular grid, grad_w_j*edge will be the central-difference,
%       so that the average (dwm+dwp)/2 will be the central-difference just
%       like in 1D. 
%--------------------------------------------------------------------------
% Edge derivative
duij = 0.5*(node(node2).u - node(node1).u);

% Left face value (wL) with the Van Albada limiter
dumL = 2*duL-duij;
dupL = duij;

% Right face value (wR) with the Van Albada limiter
dumR = -(2*duR-duij);
dupR = -duij;

switch limiter
    case 'MM'
        uL = node(node1).w + minmod([dumL,dupL]);
        uR = node(node2).w + minmod([dumR,dupR]);
    case 'VA'
        uL = node(node1).w + vanAlbada(dumL,dupL,mag_e12);
        uR = node(node2).w + vanAlbada(dumR,dupR,mag_e12);
    otherwise % Simple linear extrapolation 
        disp('no limiter set');
        uL = node(node1).w + duL;
        uR = node(node2).w - duR;
end

%  Compute the numerical flux for given uL and uR.
switch fluxMth
    case 'LF', [num_flux,wsn] = LFflux(uL,uR,n12); % Lax-Friedrichs
    case 'UP', [num_flux,wsn] = UPflux(uL,uR,n12); % Upwind flux
    otherwise 
        error('flux method not set');
end

%  Add the flux multiplied by the magnitude of the directed area vector to node1,
%  and accumulate the max wave speed quantity for use in the time step calculation.

     node(node1).res = node(node1).res + num_flux * mag_n12;
     node(node1).wsn = node(node1).wsn +      wsn * mag_n12;

%  Subtract the flux multiplied by the magnitude of the directed area vector from node2,
%  and accumulate the max wave speed quantity for use in the time step calculation.
%
%  NOTE: Subtract because the outward face normal is -n12 for the node2.

     node(node2).res = node(node2).res - num_flux * mag_n12;
     node(node2).wsn = node(node2).wsn +      wsn * mag_n12;

%-------------------------------------------------------------------------
% Close with the boundary flux using the element-based formula that is
% exact for linear fluxes (See Nishikawa AIAA2010-5093 for boundary weights
% that ensure the linear exactness for 2D/3D elements).
%
%      |  Interior Domain          |
%      |        .........          |
%      |        .       .          |
%      |        .       .          |
%      o--o--o-----o---------o--o--o  <- Boundary segment
%                  n1   |   n2
%                       v
%                     n12 (unit face normal vector)
%
% NOTE: We visit each boundary face, defined by the nodes n1 and n2,
%       and compute the flux across the boundary face: left half for node1,
%       and the right half for node2. In the above figure, the dots indicate
%       the control volume around the node n1. Clearly, the flux across the
%       left half of the face contributes to the node n1. Similarly for n2.
%
%--------------------------------------------------------------------------
%  BC: Upwind flux via freestream values
%
%      NOTE: If the final solution at the boundary node is far from
%            the freestream values, then the domain is probably is not
%            large enough. 
%--------------------------------------------------------------------------
%  BC: Solid body and Supersonic outflow
%
%      NOTE: Basically, simply compute the physical flux, which
%            can be done by calling the Roe flux with wR = wL.
%            It is equivalent to the interior-extrapolation condition.
%      NOTE: Tangency condition for solid body will be applied later.
%--------------------------------------------------------------------------
%  BC: Subsonic Outflow - Fixed Back Pressure
%
%      NOTE: Fix the pressure as freestream pressure
%            on the right side of the face (outside the domain).
%            Assumption is that the outflow boundary is far from the body.
%--------------------------------------------------------------------------

for i = 1:nbound
    for  j = 1:bound(i).nbfaces
        % Get Left and Right base data:
        n1 = bound(i).bnode( j );       %Left node
        n2 = bound(i).bnode(j+1);       %Right node
        n12(1) = bound(i).bfnx(j);      %x-component of the unit face normal vector
        n12(2) = bound(i).bfny(j);      %y-component of the unit face normal vector
        mag_e12 = bound(i).bfn(j)*0.5;  %Half length of the boundary face, j.
        switch bound(i).bcname
            case 'Neumann' % bnodes_numerical_flux_via_Neumann_BC
                %   1. Left node
                wL = node(n1).w;    wR = wL;
                [num_flux,wsn] = fluxfunc2d(wL,wR,n12,fType);
                bfluxL = num_flux;
                node(n1).wsn = node(n1).wsn + wsn*mag_e12;
                
                %   2. Right node
                wL = node(n2).w;    wR = wL;
                [num_flux,wsn] = fluxfunc2d(wL,wR,n12,fType);
                bfluxR = num_flux;
                node(n2).wsn = node(n2).wsn + wsn*mag_e12;
        end
        %   3. Add contributions to the two nodes (See Nishikawa AIAA2010-5093)
        switch elm( bound(i).belm(j) ).nV
            case 3	% Triangle
                node(n1).res = node(n1).res + (5*bfluxL+bfluxR)/6*mag_e12;
                node(n2).res = node(n2).res + (5*bfluxR+bfluxL)/6*mag_e12;
            case 4	% Quad
                node(n1).res = node(n1).res + bfluxL.mag_e12;
                node(n2).res = node(n2).res + bfluxR.mag_e12;
            otherwise
                error('Element is neither tria nor quad. Stop. ');
        end
    end
end
     
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Limiter and Numerical Fluxes  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mm = minmod(v)
    % Using Harten's generalized definition
    % minmod: zero if opposite sign, otherwise the one of smaller magnitude.
    s = sum(sign(v))/numel(v); 
    if abs(s)==1; mm = s*min(abs(v(:))); else, mm=0; end
    %m=size(v,1); mm=zeros(size(v,1),1); s=sum(sign(v),2)/m; ids=find(abs(s)==1);
    %if(~isempty(ids)); mm(ids)=s(ids).*min(abs(v(ids,:)),[],2); end
end

function va = vanAlbada(da,db,h)
    % Van Albada Slope Limiter Function
    % vanAlbada: extend the simetric formulation of the van leer limiter
    eps2=(0.3*h)^3; 
    va=0.5*(sign(da)*sign(db)+1)*((db^2+eps2)*da+(da^2+eps2)*db)/(da^2+db^2+2*eps2);
end

function [flux,swn] = LFflux(uL,uR,n12)
    % Find alpha = max|f'(u)|
    swn = max(abs(dflux(u)));

    % Compute Fj+1/2 and Fj-1/2
    for j = 2:N-1
        % Apply L-W Flux Function
        flux = 1/2*((flux(uR)+flux(uL)) - swn*(um(j+1)-up(j))); % Fj+1/2
    end
end

function [flux,swn] = UPflux(uL,uR,n12)
    % Find alpha = max|f'(u)|
    swn = max(abs(dflux(u)));

    % Compute Fj+1/2 and Fj-1/2
    for j = 2:N-1
        % Apply L-W Flux Function
        flux = 1/2*((flux(uR)+flux(uL)) - swn*(um(j+1)-up(j))); % Fj+1/2
    end
end
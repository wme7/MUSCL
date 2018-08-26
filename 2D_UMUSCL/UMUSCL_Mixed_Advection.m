%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%        Unstructured MUSCL solver for 2-d scalar advection equation
%                      by Manuel Diaz, NTU, 29.04.2015
%
%                          u_t + u_x + u_y = 0,
%
%  Node-centered finite-volume solver for the unsteady solution of the
%  linear advection equation. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Refs:
% [1] 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;

% Parameters
  CFL  = 0.25; % CFL condition;
  tEnd = 0.10; % final time;
limiter= 'VA'; % MM, VA, none;
fluxMth= 'LF'; % Lax-Friedrichs;
    nx = 11;   % number nodes in x direction;
    ny = 11;   % number nodes in y direction;

% Define a flux function
fluxfun = 'linear'; % select a flux function
switch fluxfun
    case 'linear' % linear advection
        a=-1.0; F=@(w) a*w; dF=@(w) a*ones(size(w));
        b=+1.0; G=@(w) b*w; dG=@(w) b*ones(size(w));
	case 'nonlinear' % Burgers' flux functions
        F = @(w) w.^2/2; dF = @(w) w;
        G = @(w) w.^2/2; dG = @(w) w;
end
   
% Load/build an unstructured mesh
[vx,vy,EtoV,nE,nN,BC] = SquareMesh(0,1,0,1,'MIXED',nx,ny);

% Build cells in the given unstructured mesh
[node,elem,edge,~,bound] = BuildUnstructuredMesh2d(vx,vy,EtoV,nE,nN,BC);

% Check integrity of the mesh
%CheckUnstructuredMesh2d(node,elem,edge);

% Define initial condition (IC)
IC=02;
switch IC
    case 01, f0 = @(x,y) zeros(size(x)) + 1*(x>=0.3 & x<0.7).*(y>=0.3 & y<0.7); % square jump
    case 02, f0 = @(x,y) exp(-30*((x-1/2).^2+(y-1/2).^2)); % gaussian profile
end

% Set IC in mesh
vx=[node.x]'; 
vy=[node.y]';
u0=f0(vx,vy);

% Nodes associated volume
vol = [node.vol]';

% Nodes Sum of the max wave speed multiplied by the face length
wsn = [node.wsn]';

% compute initial time step, dt0:
dt0 = ComputeTimeStep(vol,wsn);

%% Solver Loop

% Set initial time & load IC
t=0; it=0; dt=dt0; u=u0;

while t < tEnd 
    % iteration time
    if t+dt>tEnd, dt=tEnd-t; end; t=t+dt;
    
    % Runge-kutta stage 1
    L=UMUSCL_AdvRHS(u,F,dF,G,dG,node,edge,bound,limiter,fluxMth);
    us=u-dt*CFL/vol*L;
    
    % Runge-kutta stage 2
    L=UMUSCL_AdvRHS(us,F,dF,G,dG,node,edge,bound,limiter,fluxMth);
    u=0.5*(u+(us-dt*CFL/vol*L));
    
    % iteration counter
    it=it+1;
    
    % compute next step time step
    dt=ComputeTimeStep(node);
end

%% Postprocess

% Visualize initial condition
figure(1);
subplot(111); hold on;
for e = 1:nE
    patch(vx(EtoV{e}),vy(EtoV{e}),u0(EtoV{e}),mean(u0(EtoV{e}))); 
end
grid on; hold off; title(['Solution profile t=',num2str(0)]);

% Visualize solution profile at t=tEnd
subplot(112); hold on;
for e = 1:nE
    patch(vx(EtoV{e}),vy(EtoV{e}),u(EtoV{e}),mean(u(EtoV{e})));
end
grid on; hold off; title(['Solution profile t=',num2str(t)]);

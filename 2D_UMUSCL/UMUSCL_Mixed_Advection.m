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
   nEx = 10;   % number of elements in x direction;
   nEy = 10;   % number of elements in y direction;

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
[vx,vy,EtoV,nE,nN,BC] = SquareMesh(0,1,0,1,'MIXED',nEx,nEy);

% Build cells in the given unstructured mesh
[node,elem,edge,~,bound] = BuildUnstructuredMesh2d(vx,vy,EtoV,nE,nN,BC);

% Check integrity of the mesh
%CheckUnstructuredMesh2d(node,elem,edge);

% Build Cell Center's grid
x=zeros(nE,1); y=zeros(nE,1);
for i = 1:nN
    x(i)=node(i).x; y(i)=node(i).y;
end

% Build Cell Average's IC 
IC=02;
switch IC
    case 01, f0 = @(x,y) zeros(size(x)) + 1*(x>=0.3 & x<0.7).*(y>=0.3 & y<0.7); % square jump
    case 02, f0 = @(x,y) exp(-30*((x-1/2).^2+(y-1/2).^2)); % gaussian profile
end
u0=f0(x,y);

% Nodes associated volume
vol = [node.vol]';

% compute initial time step, dt0:
dt0 = ComputeTimeStep(node);

%% Solver Loop

% Set initial time & load IC
t=0; it=0; dt=dt0; u=u0;

while t < tEnd 
    % iteration time
    if t+dt>tEnd, dt=tEnd-t; end; t=t+dt;
    
    % Runge-kutta stage 1
    L=UMUSCLrhs(u,F,dF,G,dG,node,edge,bound,limiter,fluxMth);
    us=u-dt*CFL/vol*L;
    
    % Runge-kutta stage 2
    L=UMUSCLrhs(us,F,dF,G,dG,node,edge,bound,limiter,fluxMth);
    u=0.5*(u+(us-dt*CFL/vol*L));
    
    % iteration counter
    it=it+1;
    
    % compute next step time step
    dt=ComputeTimeStep(node);
end

%% Postprocess

% Visualize initial condition
figure(2); hold on;
for i = 1:nE
    patch(vx(EtoV{i})',vy(EtoV{i})',f0(vx(EtoV{i})',vy(EtoV{i})'),u0(i));
end
grid on; hold off

% Visualize solution profile
figure(3); hold on;
for i = 1:nE
    patch(vx(EtoV{i})',vy(EtoV{i})',f0(vx(EtoV{i})',vy(EtoV{i})'),u(i));
end
grid on; hold off

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%               basic MUSCL solver for scalar advection equation
%                      by Manuel Diaz, NTU, 29.04.2015
%
%                        u_t + f(u)_x + g(u)_y = 0,
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
clear; close all; clc;

for n=1:4
%% Parameters
   nE = 040;    % base number of cells/elements;
   nx = nE*2^(n-1);	% actual number of cells
   ny = nE*2^(n-1);	% actual number of cells
  CFL = 0.40;	% Courant Number;
 tEnd = 2.00;   % End time;
limiter ='MM';  % MC, MM, VA;
RKmethod='RK2'; % RK2, RK3.
plotFigs= false;

fluxfun='linear'; % select flux function
% Define our Flux function
switch fluxfun
    case 'linear'  % Scalar Advection CFL_max: 0.40
        c=1; flux = @(w) c*w; 
        dflux = @(w) c*ones(size(w));
    case 'burgers' % Burgers, CFL_max: 0.40
        flux = @(w) w.^2/2; 
        dflux = @(w) w; 
end

sourcefun='dont'; % add source term
% Source term
switch sourcefun
    case 'add', S = @(w) 0.1*w.^2;
    case 'dont',S = @(w) zeros(size(w));
end

% Build discrete domain
ax=-1; bx=1; dx=(bx-ax)/nx; xc=ax+dx/2:dx:bx;
ay=-1; by=1; dy=(by-ay)/ny; yc=ay+dy/2:dy:by;
[x,y]=meshgrid(xc,yc);

% Build IC
IC = 2;
switch IC
    case 1, u0=(x<0.5 & x>-0.5 & y<0.5 & y>-0.5); % Square block
    case 2, u0=1.0*exp(-(x.^2+y.^2)/0.1); % Smooth Guassian
end

% Plot range
plotRange=[ax,bx,ay,by,0,1];

% Initialize parpool
poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj); parpool('local',4); end

%% Solver Loop

% Load initial conditions
t=0; it=0; u=u0;

while t < tEnd
    % Update/correct time step
    dt=CFL*dx/max(abs(u(:))); if t+dt>tEnd, dt=tEnd-t; end; t=t+dt;
    
    switch RKmethod
        case 'RK2' % SSP-RK2
            % 1st stage
            L=MUSCL_AdvecRes2d_periodic(u,flux,dflux,S,dx,dy,nx,ny,limiter);	us=u-dt*L;
            
            % 2nd stage
            L=MUSCL_AdvecRes2d_periodic(us,flux,dflux,S,dx,dy,nx,ny,limiter);	u=(u+us-dt*L)/2;
            
        case 'RK3' % SSP-RK33
            % Initial step
            uo = u;

            % 1st stage
            L=MUSCL_AdvecRes2d_periodic(u,flux,dflux,S,dx,dy,nx,ny,limiter);	u=uo-dt*L;

            % 2nd Stage
            L=MUSCL_AdvecRes2d_periodic(u,flux,dflux,S,dx,dy,nx,ny,limiter);	u=0.75*uo+0.25*(u-dt*L);

            % 3rd stage
            L=MUSCL_AdvecRes2d_periodic(u,flux,dflux,S,dx,dy,nx,ny,limiter);	u=(uo+2*(u-dt*L))/3;
            
        otherwise
            error('Time integration method not set');
    end
    
    % Update iteration counter
    it=it+1;
    
    % Plot every 10 iter
    if rem(it,10)==0 && plotFigs
        surf(x,y,u,'edgecolor','none'); axis(plotRange); shg; drawnow; 
    end
    
end

%% Post Process
% Compute error norms
ue=u0; err=abs(ue(:)-u(:)); % error measurements only valid for the linear test!
L1 = dx*dy*sum(abs(err)); fprintf('L_1 norm: %1.2e \n',L1);
L2 = dx*dy*(sum(err.^2))^0.5; fprintf('L_2 norm: %1.2e \n',L2);
Linf = norm(err,inf); fprintf('L_inf norm: %1.2e \n',Linf);

% Collect measurements for convergence
Data(n).L1 = L1; %#ok<*SAGROW>
Data(n).L2 = L2;
Data(n).Linf=Linf;
Data(n).dx = dx;
Data(n).dy = dy;
Data(n).dt = dt;
if n>1 
    Data(n).rateL1Space=log(Data(n-1).L1/Data(n).L1)/log(Data(n-1).dx/Data(n).dx);
    Data(n).rateL2Space=log(Data(n-1).L2/Data(n).L2)/log(Data(n-1).dx/Data(n).dx);
end
end
disp(struct2table(Data));

%Plots results
contour(x,y,u); axis square; xlabel('x'); ylabel('y'); 
title([RKmethod,'-MUSCL Nonlinear Advection: ',fluxfun,' test 2d']);
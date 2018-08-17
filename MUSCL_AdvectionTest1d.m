%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%               basic MUSCL solver for scalar advection equation
%                      by Manuel Diaz, NTU, 29.04.2015
%
%                             u_t + f(u)_x = 0,
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
   nE = 0100;    % base number of cells;
   nx = nE*2^(n-1);	% actual number of cells
  CFL = 0.20;	 % Courant Number;
 tEnd = 2.00;    % End time;
limiter ='MC';   % MC, MM, VA;
RKmethod='RK3';  % RK2, RK3.
plotFigs= false; % plot figures

fluxfun='linear'; % select flux function
% Define our Flux function
switch fluxfun
    case 'linear'  % Scalar Advection CFL_max: 0.90
        c=1; flux = @(w) c*w; 
        dflux = @(w) c*ones(size(w));
    case 'burgers' % Burgers, CFL_max: 0.90
        flux = @(w) w.^2/2; 
        dflux = @(w) w; 
    case 'buckley' % Buckley-Leverett, CFL_max: 0.40 & tEnd: 0.40
        flux = @(w) 4*w.^2./(4*w.^2+(1-w).^2);
        dflux = @(w) 8*w.*(1-w)./(5*w.^2-2*w+1).^2;
end

sourcefun='dont'; % add source term
% Source term
switch sourcefun
    case 'add'
        S = @(w) 0.1*w.^2;
    case 'dont'
        S = @(w) zeros(size(w));
end

% Build discrete domain
a=-1; b=1; dx=(b-a)/nx; x=a+dx/2:dx:b; 

% Build IC
ICcase=2;  % {1}Testing, {2}Costum ICs
switch ICcase
    case 1 % Testing IC
        u0=TestingIC(x);  % Jiang and Shu IC
    case 2 % Guassian IC
        u0=Scalar_IC(x,7); % cases 1-9 <- check them out!
    otherwise
        error('IC file not listed');
end

% Plot range
dl=0.1; plotRange=[a,b,min(u0)-dl,max(u0)+dl];

%% Solver Loop

% Load initial conditions
t=0; it=0; u=u0;

while t < tEnd
    % Update/correct time step
    dt=CFL*dx/max(abs(u)); if t+dt>tEnd, dt=tEnd-t; end; t=t+dt;
    
    switch RKmethod
        case 'RK2' % SSP-RK2
            % 1st stage
            L=MUSCL_AdvecRes1d_periodic(u,flux,dflux,S,dx,limiter);  us=u-dt*L;
            
            % 2nd stage
            L=MUSCL_AdvecRes1d_periodic(us,flux,dflux,S,dx,limiter); u=(u+us-dt*L)/2;
            
        case 'RK3' % SSP-RK33
            % Initial step
            uo = u;

            % 1st stage
            L=MUSCL_AdvecRes1d_periodic(u,flux,dflux,S,dx,limiter);	u=uo-dt*L;

            % 2nd Stage
            L=MUSCL_AdvecRes1d_periodic(u,flux,dflux,S,dx,limiter);	u=0.75*uo+0.25*(u-dt*L);

            % 3rd stage
            L=MUSCL_AdvecRes1d_periodic(u,flux,dflux,S,dx,limiter);	u=(uo+2*(u-dt*L))/3;
            
        otherwise
            error('Time integration method not set');
    end
    
    % Update iteration counter
    it=it+1;
    
    % Plot every 10 iter
    if rem(it,10)==0 && plotFigs
        plot(x,u); axis(plotRange); shg; drawnow; 
    end
    
end

%% Post Process

% Compute error norms
ue=u0; err=abs(ue(:)-u(:)); % error measurements only valid for the linear test!
L1 = dx*sum(abs(err)); fprintf('L_1 norm: %1.2e \n',L1);
L2 = dx*(sum(err.^2))^0.5; fprintf('L_2 norm: %1.2e \n',L2);
Linf = norm(err,inf); fprintf('L_inf norm: %1.2e \n',Linf);

% Collect measurements for convergence
Data(n).L1 = L1; %#ok<*SAGROW>
Data(n).L2 = L2;
Data(n).Linf=Linf;
Data(n).dx = dx;
Data(n).dt = dt;
if n>1 
    Data(n).rateL1Space=log(Data(n-1).L1/Data(n).L1)/log(Data(n-1).dx/Data(n).dx);
    Data(n).rateL2Space=log(Data(n-1).L2/Data(n).L2)/log(Data(n-1).dx/Data(n).dx);
end
end
disp(struct2table(Data));

%Plots results
plot(x,u0,'-b',x, u,'or','MarkerSize',5); axis(plotRange);
xlabel('x'); ylabel('U'); legend('Exact', 'MUSCL'); legend boxoff;
title([RKmethod,'-MUSCL Nonlinear Advection: ',fluxfun,' test']);
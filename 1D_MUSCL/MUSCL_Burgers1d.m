%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                  basic MUSCL solver for Burger's equation
%                      by Manuel Diaz, NTU, 29.04.2015
%
%                               u_t + u*u_x = 0,
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

% Define problem constants
  A = 1.0;
cfl = 0.70;
limiter = 'MM';

% Test cases
test = 3;
switch test
    case 1
        Lx=1; nx=100; dx=Lx/nx; xc=dx/2:dx:Lx; range=[0,1,-0.25,1.25];
        U0 = 0.1*ones(size(xc)) + 1.0*(xc>=0.3 & xc<0.7); tEnd= 2;
    case 2
        Lx=2; nx=100; dx=Lx/nx; xc=dx/2:dx:Lx; range=[0,2,-0.25,1.25];
        U0 = 1.0*exp(-8*(xc-1).^2); tEnd= 2;
    case 3
        Lx=1; nx=100; dx=Lx/nx; xc=dx/2:dx:Lx; range=[0,1,-0.5,2.5];
        U0 = 1.0 + sin(2*pi*xc); tEnd= 2;
end

% Adjust grid for periodic BCs
nx=nx+2; U0=[0,U0,0];

% Boundary ghost cells
U0(1) = U0(nx-1); U0(nx) = U0(2);   % Periodic BCs

% initial time step
dt0=cfl*dx/abs(A);

% Load initial conditions
U=U0; time=0; iter=0; dt=dt0; 

% Solver Loop
while time < tEnd
    
    % Find Un+1 for 1st step
    U_s = U + dt*MUSCL_AdvecRes1d(U,dx,nx,limiter);
    
    % Find Un+2 for 2nd step
    U_s2 = U_s + dt*MUSCL_AdvecRes1d(U_s,dx,nx,limiter);
    
    % Corrector for 3rd step
    Un = 1/2*(U + U_s2);
    
    % Update U
    U = Un;
    
    % Update dt and time
    dt = cfl*dx/max(U(:));
    if time+dt>tEnd; dt=tEnd-time; end 
	time=time+dt; iter=iter+1;
    
    % Plot every 10 iter
    if rem(iter,10)==0; plot(xc,U(2:nx-1)); 
        axis(range); drawnow; end
    
end

% Remove ghost cells
U0=U0(2:nx-1); U=U(2:nx-1); nx=nx-2; 

%Plots results
plot(xc,U0,'-b',xc, U,'or','MarkerSize',5); axis(range);
xlabel('x'); ylabel('U'); legend('Exact', 'MUSCL');
title('predictor-corrector RK2 TVD-MUSCL Burgers` Eq.')
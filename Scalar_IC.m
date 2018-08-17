function u0 = Scalar_IC(x,ICcase)
% Create vector u0 with an initial condition (IC).
%**************************************************************************
% cases available: 
% {1} Gaussian wave for Advection problems
% {2} Lifted Sinusoidal wave
% {3} Centered Sinusoidal wave
% {4} HyperTangent for Burgers' equation
% {5} Riemann IC
% {6} HyperTangent
% {7} Square Jump
% {8} Displaced Square Jump for Buckley-Leverett problem
% {9} Oleg's trapezoidal
%
% Coded by Manuel Diaz 2012.12.06
%**************************************************************************
Lx=x(end)-x(1); xmid=0.5*(x(end)+x(1));
% Create the selected IC
switch ICcase
    case 1 % Gaussian wave for [-1,1] 
           % for advection test with periodic BCs. See ref [1].
        u0 = exp(-20*(x-xmid).^2);
                       
    case 2 % Centered sinusoidal wave for [-1,1] for advection.
        u0 = sin(pi*x);
        
	case 3 % Lifted sinusoidal wave for [-1,1] for advection.
        u0 = 0.5 - sin(pi*x);
        
    case 4 % Hyperbolic Tangent in [-,2,2] 
           % for Viscous Burgers' equation with Dirichlet BCs. See ref [1].
        mu = 0.02;
        u0 = 0.5*(1-tanh(x/(4*mu)));
        
    case 5 % Riemann problem
        % u = 1 for x <  x_mid
        % u = 2 for x >= x_mid
        u0 = ones(size(x));
        u0(x<=xmid) = 2;
    
    case 6 % Hyperbolic Tangent
        % u = 1 for x <  x_mid
        % u = 0 for x >= x_mid
        a = x(1); b = x(end);
        xi = (4-(-4))/(b-a)*(x - a) - 4;
        u0 = 1/2*(tanh(-4*xi)+1);
        
    case 7 % Square Jump
        u0 = rectangularPulse(xmid-0.1*Lx,xmid+0.1*Lx,x)+1;
        
	case 8 % Displaced square Jump for Buckley-Leverett problem
        xmid = -0.25;
        u0 = rectangularPulse(xmid-0.125*Lx,xmid+0.125*Lx,x);
        
    case 9 % Oleg's trapezoidal
        u0 = exp(-x).*rectangularPulse(xmid-0.1*Lx,xmid+0.1*Lx,x)*exp(.1);
        
    otherwise
        error('case not in the list')
end
function [dtU] = MUSCL_Scalar(A,U,dx,N)
%MUSCL Monotonic Upstreat Centered Scheme for Conservation Laws
%  Van Leer's MUSCL reconstruction scheme using piece wise
%  linear reconstruction
%  
%  System:
%  u_t = f(u)_x, 
%  f(u) = A*u_x,
%  f'(u) = A.
limiter = 'MC';

% Compute and limit slopes
for j = 2:N-1
    switch limiter
        case 'MC'
            % Find Sj = minmod{fwd diff, bwd diff, cntl diff} % MC limiter
            a = 2*(U(j+1) - U(j))/dx;
            b = 2*(U(j) - U(j-1))/dx;
            c = (U(j+1) - U(j-1))/(2*dx);
            S(j) = minmod([a,b,c]);
        case 'Minmod'
            % Find Sj = minmod{fwd diff, bwd diff, cntl diff} % MC limiter
            a = (U(j+1) - U(j))/dx;
            b = (U(j) - U(j-1))/dx;
            S(j) = minmod([a,b]);
    end
    
    % Build p-w cell reconstructions
    Up(j) = U(j) + S(j)*dx/2;	% Uj+1/2 -
    Um(j) = U(j) - S(j)*dx/2;	% Uj-1/2 +
    
end

% boundary conditions
Up(1) = Up(N-1); Up(N) = Up(2);   % mirrored end points
Um(1) = Um(N-1); Um(N) = Um(2);

% Find alpha = max|f'(u)|
a = max(abs(A));

% Compute Fj+1/2 and Fj-1/2
for j = 2:N-1
    
    % Apply L-W Flux Function
    Fp(j) = 1/2*((A*Up(j)+A*Um(j+1)) - a*(Um(j+1)-Up(j))); % Fj+1/2
    Fm(j) = 1/2*((A*Up(j-1)+A*Um(j)) - a*(Um(j)-Up(j-1))); % Fj-1/2
    
    % Determine d/dt Uj
    dtU(j) = -1/dx*(Fp(j) - Fm(j));
    
end

% boundary conditions
dtU(1) = dtU(N-1); dtU(N) = dtU(2);   % mirrored end points

end

function mm = minmod(v)
    % Harten's Generalized definition
    s = sum(sign(v))/numel(v); 
    if abs(s)==1; mm = s*min(abs(v(:))); else mm=0; end
end

%% DO NOT DELETE!

%     % Find Sj = minmod{fwd diff, bwd diff} 
%     a = 2*(U(j+1) - U(j))/dx;
%     b = 2*(U(j) - U(j-1))/dx;
%     S(j) = minmod(a,b);
%
% function mm = minmod(a,b)
%     % Leveque's definition
%     mm =    ((abs(a)<=abs(b)).*(a.*b>0)).*a;
%     mm = mm+((abs(b)< abs(a)).*(a.*b>0)).*b;
% end
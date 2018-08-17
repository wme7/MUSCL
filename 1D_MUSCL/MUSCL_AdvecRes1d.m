function [res] = MUSCL_AdvecRes1d(u,dx,N,limiter)
%MUSCL Monotonic Upstreat Centered Scheme for Conservation Laws
%  Van Leer's MUSCL reconstruction scheme using piece wise
%  linear reconstruction

% for burgers equation
flux = @(w) w.^2/2; 
dflux = @(w) w;

% Compute and limit slopes
du=zeros(1,N); up=zeros(1,N); um=zeros(1,N);
for j = 2:N-1
    switch limiter
        case 'MC'
            % Find Sj = minmod{fwd diff, bwd diff, cntl diff} % MC limiter
            a = 2*(u(j+1) - u(j))/dx;
            b = 2*(u(j) - u(j-1))/dx;
            c = (u(j+1) - u(j-1))/(2*dx);
            du(j) = minmod([a,b,c]);
        case 'MM'
            % Find Sj = minmod{fwd diff, bwd diff, cntl diff} % MC limiter
            a = (u(j+1) - u(j))/dx;
            b = (u(j) - u(j-1))/dx;
            du(j) = minmod([a,b]);
        case 'VA'
            a = (u(j+1) - u(j))/dx;
            b = (u(j) - u(j-1))/dx;
            du(j) = vanalbada(a,b,dx);
    end
    
    % Build p-w cell reconstructions
    up(j) = u(j) + du(j)*dx/2;	% U_j+1/2^{-}
    um(j) = u(j) - du(j)*dx/2;	% U_j-1/2^{+}
    
end

% boundary conditions
up(1) = up(N-1); up(N) = up(2);   % mirrored end points
um(1) = um(N-1); um(N) = um(2);

% Find alpha = max|f'(u)|
a = max(abs(dflux(u)));

% Compute Fj+1/2 and Fj-1/2
for j = 2:N-1
    
    % Apply L-W Flux Function
    Fp(j) = 1/2*((flux(up(j))+flux(um(j+1))) - a*(um(j+1)-up(j))); % Fj+1/2
    Fm(j) = 1/2*((flux(up(j-1))+flux(um(j))) - a*(um(j)-up(j-1))); % Fj-1/2
    
    % Determine d/dt Uj
    res(j) = -(Fp(j)-Fm(j))/dx;
    
end

% boundary conditions
res(1) = res(N-1); res(N) = res(2);   % mirrored end points

end

function mm = minmod(v)
    % Using Harten's generalized definition
    % minmod: zero if opposite sign, otherwise the one of smaller magnitude.
    %m=size(v,1); mm=zeros(size(v,2),1); s=sum(sign(v),2)/m; ids=find(abs(s)==1);
    %if(~isempty(ids)); mm(ids)=s(ids).*min(abs(v(ids,:)),[],2); end
    s = sum(sign(v))/numel(v); 
    if abs(s)==1; mm = s*min(abs(v(:))); else, mm=0; end
end

function va = vanalbada(da,db,h)
    % Van Albada Slope Limiter Function
    % vanAlbada: extend the simetric formulation of the van leer limiter
    eps2=(0.3*h)^3; 
    va=0.5*(sign(da)*sign(db)+1)*((db^2+eps2)*da+(da^2+eps2)*db)/(da^2+db^2+2*eps2);
end

%% DO NOT DELETE!
% Some references I used in my development.
%
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
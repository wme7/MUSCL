function [res] = MUSCL_AdvecRes2d_periodic(qi,flux,dflux,S,dx,dy,N,M,limiter)
%   MUSCL Monotonic Upstreat Centered Scheme for Conservation Laws
%   Van Leer's MUSCL reconstruction scheme using piece wise linear
%   reconstruction
%
% Written by Manuel Diaz, NTU, 04.29.2015.

% Normal unitary face vectors: (nx,ny)
%normals = {[0,1], [1,0], [0,-1], [-1,0]}; % i.e.: [N, E, S, W]

% Initial Arrays      
% qi : q_{ i,j }^{n},
qim1 = circshift(qi,[0, 1]); % : q_{j,i-1}^{n},
qip1 = circshift(qi,[0,-1]); % : q_{j,i+1}^{n}.
qjm1 = circshift(qi,[ 1,0]); % : q_{j-1,i}^{n}.
qjp1 = circshift(qi,[-1,0]); % : q_{j+1,i}^{n}.

dqE = qip1-qi; dqW = qi-qim1; dqC = (qip1-qim1)/2; dqdx = zeros(size(qi));
dqN = qjp1-qi; dqS = qi-qjm1; dqM = (qjp1-qjm1)/2; dqdy = zeros(size(qi)); 

% Compute and limit slopes at cells (i,j)
parfor i = 1:M % only internal cells
    for j = 1:N % only internal cells
        switch limiter
            case 'MC' % MC limiter: Find dq_j = minmod{fwd diff, bwd diff, cntrl diff}
                dqdx(j,i) = minmod([dqW(j,i),dqE(j,i),dqC(j,i)]);
                dqdy(j,i) = minmod([dqS(j,i),dqN(j,i),dqM(j,i)]);
            case 'MM' % Minmod limiter: Find dq_j = minmod{fwd diff, bwd diff}
                dqdx(j,i) = minmod([dqW(j,i),dqE(j,i)]);
                dqdy(j,i) = minmod([dqS(j,i),dqN(j,i)]);
            case 'VA' % Van Albada limiter
                dqdx(j,i) = vanalbada(dqW(j,i),dqE(j,i),dx);
                dqdy(j,i) = vanalbada(dqS(j,i),dqN(j,i),dy);
        end
    end
end

%% Compute MUSCL reconstructions

% Left and Right extrapolated q-values at the boundary i+1/2
qiph = qi+dqdx/2;	% q_{j,i+1/2}^{-} from cell_ij
qimh = qi-dqdx/2;	% q_{j,i-1/2}^{+} from cell_ij

qL=qiph; qR=circshift(qimh,[0,-1]);

% Compute Lax-Friedrichs numerical flux and update solution
LFx = 0.5*(flux(qL)+flux(qR)-abs(dflux((qi+qip1)/2)).*(qR-qL)); % Lax friedrichs flux

qiph = qi+dqdy/2;	% q_{j+1/2,i}^{-} from cell_ij
qimh = qi-dqdy/2;	% q_{j-1/2,i}^{+} from cell_ij

qL=qiph; qR=circshift(qimh,[-1,0]);

% Compute Lax-Friedrichs numerical flux and 
LFy = 0.5*(flux(qL)+flux(qR)-abs(dflux((qi+qjp1)/2)).*(qR-qL)); % Lax friedrichs flux

% Compute operator L = - df(q)/dx - dg(q)/dy + S(q).
res = (LFx-circshift(LFx,[0,1]))/dx + (LFy-circshift(LFy,[1,0]))/dy - S(qi); % 

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
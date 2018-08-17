function [res] = MUSCL_AdvecRes1d_periodic(qi,flux,dflux,S,dx,limiter)
%MUSCL Monotonic Upstreat Centered Scheme for Conservation Laws
%  Van Leer's MUSCL reconstruction scheme using piece wise
%  linear reconstruction
%
% Written by Manuel Diaz, NTU, 04.29.2015.

% Initial Arrays      
% qi = q;  % : q_{ j }^{n},
qim1 = circshift(qi,+1); % : q_{j-1}^{n},
qip1 = circshift(qi,-1); % : q_{j+1}^{n}.
 
dqR = qip1 - qi;
dqL = qi - qim1;
dqC = (qip1-qim1)/2;
dq=zeros(size(qi));

% Compute and limit slopes
for j = 1:size(qi,2) % for all internal faces
    switch limiter
        case 'MC',dq(j) = minmod([2*dqR(j),2*dqL(j),dqC(j)]); % MC limiter: Find dq_j = minmod{fwd diff, bwd diff, cntrl diff}
        case 'MM',dq(j) = minmod([dqR(j),dqL(j)]); % Minmod limiter: Find dq_j = minmod{fwd diff, bwd diff}
        case 'VA',dq(j) = vanalbada(dqR(j),dqL(j),dx); % Van Albada limiter.
        otherwise, error('limiter not listed!');
    end
end

% Left and Right extrapolated q-values at the boundary j+1/2
qiph_M = qi+dq/2;	% q_{j+1/2}^{-} from cell j
qimh_M = qi-dq/2;	% q_{j+1/2}^{+} from cell j

qL=circshift(qiph_M,0);
qR=circshift(qimh_M,-1);

%% Compute Lax-Friedrichs numerical flux and update solution
LF = 0.5*(flux(qL)+flux(qR)-abs(dflux((qi+qip1)/2)).*(qR-qL)); % Lax friedrichs flux
res = (LF-circshift(LF,1))/dx - S(qi); % L = df(q)/dx + S(q).

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
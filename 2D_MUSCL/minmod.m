function mm = minmod(v)
    % Harten's Generalized definition
    s = sum(sign(v))/numel(v); 
    if abs(s)==1; mm = s*min(abs(v(:))); else, mm=0; end
end
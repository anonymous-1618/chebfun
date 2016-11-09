function rts = findExtrema(f, p, q, rh, h, xk)
% finds all the local maxima and minima
% of f-p./q
% xk is the present reference
% rh is a handle to p/q
% h is the current leveled error

err_handle = @(x) feval(f, x) - rh(x);
rts = [];

doms = unique([f.domain(1); xk; f.domain(end)]);

% Initial trial
if ( isempty(xk) )
    
    % TODO: there is generally a warning here that the function
    % was not resolved when dealing with more complicated instances
    % (i.e., splitting on was required to successfuly work with f).
    % Can we eliminate it?
    ek = chebfun(@(x) err_handle(x), doms.', 'splitting', 'on');
    rts = roots(diff(ek), 'nobreaks');

    % The following lines failed sometimes to find all the local extrema
    % of the error function of f, when f is composed of several pieces
    
    %e_num = (q.^2).*diff(f) - q.*diff(p) + p.*diff(q);
    %rts = roots(e_num, 'nobreaks');

    
    rts = unique([f.domain(1); rts; f.domain(end)]);
end
   
if ( ~isempty(xk) )
    for k = 1:length(doms)-1
        ek = chebfun(@(x) err_handle(x), [doms(k), doms(k+1)], 24, 'eps', 1e-12); 
        %plotcoeffs(ek);
        %pause
        ek = simplify(ek);
        rts = [rts; roots(diff(ek), 'nobreaks')];  %#ok<AGROW>
    end    
end

% Append end points of the domain:
rts = [f.domain(1) ; rts; f.domain(end)];

end


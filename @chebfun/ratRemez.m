function varargout = ratRemez(varargin)
%REMEZ   Best polynomial or rational approximation for real valued chebfuns.
%   P = REMEZ(F, M) computes the best polynomial approximation of degree M to
%   the real CHEBFUN F in the infinity norm using the Remez algorithm.
%
%   [P, Q] = REMEZ(F, M, N) computes the best rational approximation P/Q of type
%   (M, N) to the real CHEBFUN F using the Remez algorithm.
%
%   [P, Q, R_HANDLE] = REMEZ(F, M, N) does the same but additionally returns a
%   function handle R_HANDLE for evaluating the rational function P/Q.
%
%   [...] = REMEZ(..., 'tol', TOL) uses the value TOL as the termination
%   tolerance on the increase of the levelled error.
%
%   [...] = REMEZ(..., 'display', 'iter') displays output at each iteration.
%
%   [...] = REMEZ(..., 'maxiter', MAXITER) sets the maximum number of allowable
%   iterations to MAXITER.
%
%   [...] = REMEZ(..., 'plotfcns', 'error') plots the error after each iteration
%   while the algorithm executes.
%
%   [P, ERR] = REMEZ(...) and [P, Q, R_HANDLE, ERR] = REMEZ(...) also returns
%   the maximum error ERR.
%
%   [P, ERR, STATUS] = REMEZ(...) and [P, Q, R_HANDLE, ERR, STATUS] = REMEZ(...)
%   also return a structure array STATUS with the following fields:
%      STATUS.DELTA  - Obtained tolerance.
%      STATUS.ITER   - Number of iterations performed.
%      STATUS.DIFFX  - Maximum correction in last trial reference.
%      STATUS.XK     - Last trial reference on which the error equioscillates.
%
%   This code is quite reliable for polynomial approximations but rather
%   fragile for rational approximations.  Better results can often be obtained
%   with CF(), especially if f is smooth.
%
% References:
%
%   [1] Pachon, R. and Trefethen, L. N.  "Barycentric-Remez algorithms for best
%   polynomial approximation in the chebfun system", BIT Numerical Mathematics,
%   49:721-742, 2009.
%
%   [2] Pachon, R.  "Algorithms for Polynomial and Rational Approximation".
%   D. Phil. Thesis, University of Oxford, 2010 (Chapter 6).
%
% See also CF.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin < 2 ) 
    error('CHEBFUN:CHEBFUN:remez:nargin', ...
        'Not enough input arguments.');
else
    [f, m, n, opts] = parseInput(varargin{:});
end

if ( isempty(f) ) 
    varargout = {f};
    return;
end

if ( opts.rationalMode )
    [m, n] = adjustDegreesForSymmetries(f, m, n);
end

dom = f.domain([1, end]);
normf = norm(f);

% With zero denominator degree, the denominator polynomial is trivial.
if ( n <= 0 ) 
    error('CHEBFUN:CHEBFUN:ratRemez:ratOnly', ...
        'Use remez for polynomial case' );
end

% Initial values for some parameters.
iter = 0;                 % Iteration count.
deltaLevelError = max(normf, eps);  % Value for stopping criterion.
deltamin = inf;           % Minimum error encountered.
deltaReference = 1;                % Maximum correction to trial reference.

N = m + n;
% Compute an initial reference set to start the algorithm.
%[xk, pmin, qmin] = getInitialReference(f, m, n, N);

% abs(x) for (12,12)
% xk = [-9.095373966236864022202741740653e-01;-6.985246094149314640205350528505e-01;-4.740095023485303225357337630444e-01;-2.964071840593504911878876963069e-01;-1.751337137726280472704422780010e-01;-9.869274440969609892819818200507e-02;-5.297163444127219010018353388894e-02;-2.683249941873676575145180088395e-02;-1.260393508312235295490110089999e-02;-5.319357635375943506999272802597e-03;-1.892241268243602084746504296364e-03;-4.794997616838592125380586569721e-04;3.500513251699756604892707340536e-18;4.794997616838592448102587123069e-04;1.892241268243601829183712777340e-03;5.319357635375943690750648083533e-03;1.260393508312235325710079105840e-02;2.683249941873676587454119586356e-02;5.297163444127219088707087985204e-02;9.869274440969610158965260038120e-02;1.751337137726280516593253572205e-01;2.964071840593504913722299166664e-01;4.740095023485303232483560172590e-01;6.985246094149314624074972612828e-01;9.095373966236864011860941942525e-01;1.000000000000000000000000000000e+00;]

% abs(x) for (14,14)
%xk = [-9.211540154753614829520378117407e-01;-7.317375192530093183445880750584e-01;-5.201393332354511588265208544860e-01;-3.431398664861511955585953603370e-01;-2.152489519894858548860885161167e-01;-1.299014584126534657955344267697e-01;-7.564540785160607506290606738029e-02;-4.239477594634671261729444303458e-02;-2.271082156244639250291487899161e-02;-1.149685102863935825219296754920e-02;-5.399378547111094898010284219310e-03;-2.278640873188680588326146261605e-03;-8.105663955549065768266619896281e-04;-2.053996505412981690390701901139e-04;-1.799141134250241541594846236570e-18;2.053996505412981450454282000124e-04;8.105663955549065519134007239022e-04;2.278640873188680568430267056384e-03;5.399378547111094900792290263712e-03;1.149685102863935823182784732886e-02;2.271082156244639254769100517534e-02;4.239477594634671259752142524997e-02;7.564540785160607543420046289672e-02;1.299014584126534656304714626971e-01;2.152489519894858531233831668601e-01;3.431398664861511935253272744408e-01;5.201393332354511609606246927205e-01;7.317375192530093105803980311560e-01;9.211540154753614896249548513299e-01;1.000000000000000000000000000000e+00];

%xk = [-1.000000000000000000000000000000e+00;-9.591701232116413120115671732307e-01;-8.506316071677443759553102923285e-01;-7.064081519061093562991422754420e-01;-5.579525507856667902112575707894e-01;-4.251052139991369730841265238238e-01;-3.157669985609909769217128871396e-01;-2.302773741582663299344531488179e-01;-1.655583999867939564766371513731e-01;-1.176006486828651114171103656622e-01;-8.260711822897671375073397608934e-02;-5.738786880449465308086150500677e-02;-3.941202685080791473882854202028e-02;-2.673187704728606289918258717435e-02;-1.787873270734587332797134650992e-02;-1.176021878511486850444855671538e-02;-7.583838830232069504167959887711e-03;-4.805189123766551606206014104730e-03;-3.017590354978343821419753668125e-03;-1.863150940469103132308034954325e-03;-1.090193760462791825007094727274e-03;-6.085253323058575768414933987577e-04;-3.386586725530730543476906970755e-04;-1.845495103313738111076323138382e-04;-9.348592763181898475491519206970e-05;-4.368568137818709755648263330349e-05;-1.817910130813921588648596965571e-05;-6.261934857463058514930242066636e-06;-1.544719411595569960182358729392e-06;-3.076031556019783504880584144403e-19;1.544719411595566483265382661199e-06;6.261934857463042872538445491475e-06;1.817910130813916843210153475218e-05;4.368568137818689069777206026864e-05;9.348592763181878798409036964426e-05;1.845495103313744437513959689272e-04;3.386586725530735384438880002861e-04;6.085253323058581045278020470072e-04;1.090193760462796786740370833536e-03;1.863150940469108173288096293769e-03;3.017590354978333836498699850109e-03;4.805189123766522524060815953795e-03;7.583838830232003585317497634473e-03;1.176021878511474432425673579626e-02;1.787873270734574335288448226174e-02;2.673187704728607943480303565494e-02;3.941202685080822466127730372034e-02;5.738786880449526449137492461893e-02;8.260711822897746983666526409033e-02;1.176006486828658165094963138211e-01;1.655583999867945061986817660436e-01;2.302773741582667192345550075499e-01;3.157669985609912403853322568231e-01;4.251052139991371465202889691344e-01;5.579525507856669026037173340776e-01;7.064081519061094295359981905863e-01;8.506316071677444196764417667185e-01;9.591701232116413485763994335364e-01;    ]



% abs(x) for (28,28)
%xk = [-1.000000000000000000000000000000e+00;-9.584370705751226040221312908452e-01;-8.481926267986693568740057512253e-01;-7.022708052382531942288155292692e-01;-5.527903770798716595562907303047e-01;-4.197069313562121550054809905074e-01;-3.107169142573433071730768774372e-01;-2.258945144540723380802586285839e-01;-1.619551490168538807505297657804e-01;-1.147611399608775188976927370737e-01;-8.044845604774520775089622634791e-02;-5.579968405119176023665255381152e-02;-3.828129835489647181601657041132e-02;-2.595863009899677569865152640713e-02;-1.738192897331838068533503577038e-02;-1.147906142869793465195863973906e-02;-7.465527594442907813603337821413e-03;-4.772815178224992867284973538717e-03;-2.992867765464625972784829949574e-03;-1.835724442718635580125258367109e-03;-1.097562341299844895072513177818e-03;-6.368042510032572383455369303557e-04;-3.564057793248883647702084701610e-04;-1.908355252374004781488550963665e-04;-9.659144031230126990718210335365e-05;-4.536113080296000798773281923190e-05;-1.914303760745804757301344991167e-05;-6.809613064518808784350406010121e-06;-1.725573055193166659042118946800e-06;-3.076031556019783504880584144403e-19;1.725573055193166460953849619666e-06;6.809613064518809086740675081301e-06;1.914303760745804681436261765488e-05;4.536113080296001028782041119461e-05;9.659144031230126539707875396038e-05;1.908355252374004710251900609309e-04;3.564057793248883751294181401172e-04;6.368042510032572541092625871539e-04;1.097562341299844834159875594932e-03;1.835724442718635613101283579415e-03;2.992867765464626003782850330997e-03;4.772815178224992852902537065986e-03;7.465527594442908010333222843773e-03;1.147906142869793466673621198754e-02;1.738192897331838060306250769630e-02;2.595863009899677555288720350763e-02;3.828129835489647191280811401656e-02;5.579968405119176002720611817117e-02;8.044845604774520761694014342178e-02;1.147611399608775192391424021117e-01;1.619551490168538807193190900543e-01;2.258945144540723355737360542104e-01;3.107169142573433082855942266888e-01;4.197069313562121561450085352107e-01;5.527903770798716557599243479781e-01;7.022708052382531919927958462886e-01;8.481926267986693310450573877954e-01;9.584370705751226054957525429332e-01];

% abs(x) for (30,30)
%xk = [-9.610645113969592789377426846687e-01;-8.570333473840675335201706853515e-01;-7.175420080195055688557359500836e-01;-5.723184844495258137836444181154e-01;-4.407710494500139112934854437877e-01;-3.311811850452570494546287776679e-01;-2.444777796824113557363126164766e-01;-1.780757946384837450195879323874e-01;-1.282943185009028267920729515933e-01;-9.152806307233569025273426479225e-02;-6.468654703344661746732194520171e-02;-4.528363776093735665393937193625e-02;-3.138701408458068849156840948042e-02;-2.152540019837480240581131461442e-02;-1.459386954686970286832932035522e-02;-9.771246057732454819272105515346e-03;-6.452692771916802445456474211179e-03;-4.196498670950387596047908532756e-03;-2.682857428789305947773421574951e-03;-1.682321251189870669557151196799e-03;-1.031877702396002058188513761225e-03;-6.169495681466759574949031612091e-04;-3.579532549150287423333490499161e-04;-2.003387966750464176269941267768e-04;-1.072703096203667702655424359322e-04;-5.429488872025163656217291696610e-05;-2.549788609111936897967099984127e-05;-1.076046790212313075212167219418e-05;-3.827742717456592090320219358674e-06;-9.699596192353484331145364556938e-07;-1.197571250825071028973076045888e-23;9.699596192353484021277893074039e-07;3.827742717456592152212060050174e-06;1.076046790212312988483230233192e-05;2.549788609111936878839429158019e-05;5.429488872025163718748607388354e-05;1.072703096203667693209433811937e-04;2.003387966750464150160854888382e-04;3.579532549150287387036211484956e-04;6.169495681466759783687781361355e-04;1.031877702396002052129500372459e-03;1.682321251189870588985862690367e-03;2.682857428789305935035833107296e-03;4.196498670950387561623041478295e-03;6.452692771916802397737342641674e-03;9.771246057732454837804695738507e-03;1.459386954686970293113521014797e-02;2.152540019837480241847311096173e-02;3.138701408458068851953396617534e-02;4.528363776093735613848576170335e-02;6.468654703344661685191974800188e-02;9.152806307233569055076370287350e-02;1.282943185009028264207045339783e-01;1.780757946384837439178568496931e-01;2.444777796824113545779005017742e-01;3.311811850452570464571720513335e-01;4.407710494500139101144281996017e-01;5.723184844495258145187422681724e-01;7.175420080195055731339341956137e-01;8.570333473840675278620480561637e-01;9.610645113969592813126916158988e-01;1.000000000000000000000000000000e+00];

% abs(x) for higher degree
xk = [logspace(-7,0,m+1)]; xk = [xk(1:end-1) 0 -xk]; xk = sort(xk,'ascend'); xk = xk';
%xk = chebpts(m+n+2);


%{
if ( isempty(qmin) )
    qmin = chebfun(1, f.domain([1, end]));
end
%}


if exist('xk','var')==0 % get initial guess from polyremez
[p,~,status] = remez(f,m+n);xk = status.xk;
end


% Print header for text output display if requested.
if ( opts.displayIter )
    disp('It.     Max(|Error|)       |ErrorRef|      Delta ErrorRef      Delta Ref')
end


% Old reference and levelled error
x0 = xk;
h0 = 0;
h = -1; hpre = 0;
% Run the main algorithm.
%while ( (deltaLevelError/normf > opts.tol) && (iter < opts.maxIter) && (deltaReference > 0) )
while ( (abs(abs(h)-abs(hpre))/abs(h) > opts.tol) && (iter < opts.maxIter) && (deltaReference > 0) )
    %disp(abs(h)-abs(hpre))    
    hpre = h; 
    [p, q, rh, pqh, h, interpSuccess, xsupport] = computeTrialFunctionRational(f, xk, m, n);      
    
       xx = linspace(-1,1,5000);
       chebfunsetting
       subplot(2,1,1)

       fxx = feval(f,xx);
       fxk = feval(f,xk);
        hold off,plot(xx,fxx-rh(xx),LW,lw),hold on
        plot(xk,fxk-rh(xk),'k.',MS,ms),drawnow
        
    subplot(2,1,2)
    xx = logspace(-8,0,1000);
    hold off
    fxx = feval(f,xx);
    semilogx(xx,fxx-rh(xx),LW,lw),hold on
    semilogx(xk(find(xk>0)),fxk(find(xk>0))-rh(xk(find(xk>0))),'k.',MS,ms),drawnow   
        
        
     
    %if ( abs(h) <= abs(h0) )
        % The levelled error has not increased
        %disp('level error decreased' )
        %xk = makeNewReference(xkPrev, xk);
    %end
    
    %perform modiffications to the previous reference set until the current
    %interpolant has no poles on the approximation domain (technically,
    %we're not doing that just yet!!!)
    if (interpSuccess == 0)
        disp('perturb current reference set');
        % xk is the current reference set that is causing us problems
        [~, pos] = max(abs(feval(errPrev, xk)));
        pos = pos(1);
        [~, closeIdx] = min(abs(xkPrev-xk(pos)));
        if (sign(errPrev(xk(pos))) ~= sign(errPrev(xkPrev(closeIdx))))
            if(closeIdx == 1)
                closeIdx = 2;
            elseif(closeIdx == length(xkPrev))
                closeIdx = length(xkPrev) - 1;
            else
                if(abs(xkPrev(closeIdx-1)-xk(pos)) < abs(xkPrev(closeIdx+1)-xk(pos)))
                    closeIdx = closeIdx - 1;
                else
                    closeIdx = closeIdx + 1;
                end    
            end
        end
        
        xkNew = xkPrev;
        xkNew(closeIdx) = xk(pos);
        [p, q, rh, pqh, h, interpSuccess] = computeTrialFunctionRational(f, xkNew, m, n);
        
        
        modifCounter = 0;
        while (interpSuccess == 0 & (modifCounter <12))
            xkNew = xkPrev;
            modifCounter = modifCounter+1;
            xkNew(closeIdx) = (1/2^modifCounter)*xk(pos) + (1 - 1/2^modifCounter)*xkPrev(closeIdx);
            [p, q, rh, pqh, h, interpSuccess] = computeTrialFunctionRational(f, xkNew, m, n);
        end
        if(interpSuccess == 0)
            error('Failed in the perturbation of the reference set!');
        end
        % Again, optimistic assumption that the previous loop works
        xk = xkNew;
    end
    
    if ( abs(h) <= abs(h0) )
        % The levelled error has not increased
        disp('level error decreased' )
        %xk = makeNewReference(xkPrev, xk);
    end
    
    xkPrev = xk;
    %hPrev = h;
    
    % Perturb exactly-zero values of the levelled error.
    if ( h == 0 )
        h = 1e-19;
    end
    %xkPrev = xk;
    % Update the exchange set using the Remez algorithm with full exchange.   
    [xk, err, err_handle] = exchange(xk, h, 2, f, p, q, rh, N + 2, opts);
    errPrev = err_handle;

    % Update max. correction to trial reference and stopping criterion.
    deltaReference = max(abs(x0 - xk));
    deltaLevelError = err - abs(h);

    % Store approximation with minimum norm.
    if ( deltaLevelError < deltamin )
        pmin = p;
        if ( n > 0 )
            qmin = q;
        end

        errmin = err;
        xkmin = xk;
        deltamin = deltaLevelError;
    end

    % Display diagnostic information as requested.
    if ( opts.plotIter )
        doPlotIter(x0, xk, err_handle, dom);
        pause();
    end

    if ( opts.displayIter )
        doDisplayIter(iter, err, h, deltaLevelError, normf, deltaReference);
    end
    
    if ( opts.demoMode )
        pause
    end

    % Save the old reference and
    % the old levelled error.
    x0 = xk;
    h0 = h;
    
    iter = iter + 1;
end

%{
% Take best results of all the iterations we ran.
p = pmin; 
err = errmin;
xk = xkmin;
deltaLevelError = deltamin;
%}

% Warn the user if we failed to converge.
if ( deltaLevelError/normf > opts.tol )
    warning('CHEBFUN:CHEBFUN:remez:convergence', ...
        ['Remez algorithm did not converge after ', num2str(iter), ...
         ' iterations to the tolerance ', num2str(opts.tol), '.']);
end

% Form the outputs.
status.delta = deltaLevelError/normf;
status.iter = iter;
status.deltaReference = deltaReference;
status.xk = xk;
status.xsupport = xsupport;

%{
p = simplify(p);
if ( opts.rationalMode )
    q = simplify(qmin);
    varargout = {p, q, @(x) feval(p, x)./feval(q, x), err, status};
else
    varargout = {p, err, status};
end
%}

if nargout<=2
    varargout = {rh,status};
else
    varargout = {p,q,rh,status};    
end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions implementing the core part of the algorithm.

function [xk, p, q] = getInitialReference(f, m, n, N)

% If doing rational Remez, get initial reference from trial function generated
% by CF or Chebyshev-Pade.
flag = 0;
a = f.domain(1);
b = f.domain(end);

if ( numel(f.funs) == 1 )
    %[p, q] = chebpade(f, m, n);
    [p, q] = cf(f, m, n);
else
    %[p, q] = chebpade(f, m, n, 5*N);
    [p, q] = cf(f, m, n, 100*N);
end
pqh = @(x) feval(p, x)./feval(q, x);
[xk, err, e, flag] = exchange([], 0, 2, f, p, q, pqh, N + 2);

% If the above procedure failed to produce a reference
% with enough equioscillation points, just use the Chebyshev points.
if ( flag == 0 )
    disp('CF failed');
    xk = chebpts(N + 2, f.domain([1, end]), 1);
    xk = [xk(round(length(xk)/2)+1:length(xk)) - 1; xk(1:round(length(xk)/2))+1];
    xk(round(length(xk)/2))=-1;
    xk = sort(xk);    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions for displaying diagnostic information.

% Function called when opts.plotIter is set.
function doPlotIter(xo, xk, err_handle, dom)

xxk = linspace(dom(1), dom(end), max(3000, 10*length(xk)));
plot(xo, 0*xo, '.r', 'MarkerSize', 6)   % Old reference.
holdState = ishold;
hold on
plot(xk, 0*xk, 'ok', 'MarkerSize', 3)   % New reference.

plot(xxk, err_handle(xxk))               % Error function.
if ( ~holdState )                        % Return to previous hold state.
    hold off
end
xlim(dom)
legend('Current Ref.', 'Next Ref.', 'Error')
drawnow

end

% Function called when opts.displayIter is set.
function doDisplayIter(iter, err, h, delta, normf, deltaReference)

disp([num2str(iter), '        ', num2str(err, '%5.4e'), '        ', ...
    num2str(abs(h), '%5.4e'), '        ', ...
    num2str(delta/normf, '%5.4e'), '        ', num2str(deltaReference, '%5.4e')])
end


function status = remezParseFunction(f)
% Parse a Chebfun to see if Remez
% can be applied to it

% Add conditions needed on f:
if ( ~isreal(f) )
    error('CHEBFUN:CHEBFUN:remez:real', ...
        'REMEZ only supports real valued functions.');    
end

if ( numColumns(f) > 1 )
    error('CHEBFUN:CHEBFUN:remez:quasi', ...
        'REMEZ does not currently support quasimatrices.');    
end

if ( issing(f) )
    error('CHEBFUN:CHEBFUN:remez:singularFunction', ...
        'REMEZ does not currently support functions with singularities.');
end

% If all are satisifed, we can go ahead:
status = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input parsing.

function [f, m, n, opts] = parseInput(varargin)

%% Handle the first two arguments:
f = varargin{1};
remezParseFunction(f);
varargin(1) = [];


m = varargin{1};
varargin(1) = [];

n = 0;
opts.rationalMode = false;
if ( ~isempty(varargin) > 0 )
    if ( isa(varargin{1}, 'double') )
        n = varargin{1};
        varargin(1) = [];
        opts.rationalMode = true;
    end
end

if ( isa(m, 'vector') || isa(n, 'vector') || ...
        m < 0 || n < 0 || m ~= round(m) || n ~= round(n) )
   error('CHEBFUN:CHEBFUN:remez:degree', ...
        'Approximant degree must be a non-negative integer.');
end


% Parse name-value option pairs.
N = m + n;
%opts.tol = 1e-16*(N^2 + 10); % Relative tolerance for deciding convergence.
opts.tol = 1e-5;
opts.maxIter = 10+round(max(m,n)/2);           % Maximum number of allowable iterations.
opts.displayIter = false;    % Print output after each iteration.
opts.plotIter = false;       % Plot approximation at each iteration.
opts.demoMode = false;

while ( ~isempty(varargin) )
    if ( strcmpi('tol', varargin{1} ) )
        opts.tol = varargin{2};
        varargin(1) =[];
        varargin(1) =[];        
    elseif ( strcmpi('maxiter', varargin{1}) )
        opts.maxIter = varargin{2};        
        varargin(1) =[];
        varargin(1) =[];        
    elseif ( strcmpi('display', varargin{1}) )
        varargin(1) = [];
        if ( strcmpi('iter', varargin{1}) )
            varargin(1) = [];
        end
        opts.displayIter = true;        
    elseif ( strcmpi('demo', varargin{1}) )        
        varargin(1) = [];
        opts.demoMode = true;
        opts.displayIter = true;
        opts.plotIter = true;        
    elseif ( strcmpi('plotfcns', varargin{1}) )
        varargin(1) = [];
        if ( strcmpi('error', varargin{1}) )
            varargin(1) = [];
        end            
        opts.plotIter = true;        
    else
        error('CHEBFUN:CHEBFUN:remez:badInput', ...
            'Unrecognized sequence of input parameters.')
    end
        
end

end


function [m, n] = adjustDegreesForSymmetries(f, m, n)
%ADJUSTDEGREESFORSYMMETRIES   Adjust rational approximation degrees to account
%   for function symmetries.
%
%   [M, N] = ADJUSTDEGREESFORSYMMETRIES(F, M, N) returns new degrees M and N to
%   correct the defect of the rational approximation if the target function is
%   even or odd.  In either case, the Walsh table is covered with blocks of
%   size 2x2, e.g.  for even function the best rational approximant is the same
%   for types [m/n], [m+1/n], [m/n+1] and [m+1/n+1], with m and n even. This
%   strategy is similar to the one proposed by van Deun and Trefethen for CF
%   approximation in Chebfun (see @chebfun/cf.m).

% Sample piecewise-smooth CHEBFUNs.
if ( (numel(f.funs) > 1) || (length(f) > 128) )
  f = chebfun(f, f.domain([1, end]), 128);
end

% Compute the Chebyshev coefficients.
c = chebcoeffs(f, length(f));
c(1) = 2*c(1);

% Check for symmetries and reduce degrees accordingly.
if ( max(abs(c(2:2:end)))/vscale(f) < eps )   % f is even.
    if ( mod(m, 2) == 1 )
        m = max(m - 1, 0);
    end
    if ( mod(n, 2) == 1 )
        n = max(n - 1, 0);
    end
elseif ( max(abs(c(1:2:end)))/vscale(f) < eps ) % f is odd.
    if ( mod(m, 2) == 0 )
        m = max(m - 1, 0);
    end
    if ( mod(n, 2) == 0 )
        n = max(n + 1, 0);
    end
end

end


function [xk, norme, err_handle, flag] = exchange(xk, h, method, f, p, q, rh, Npts, opts)
%EXCHANGE   Modify an equioscillation reference using the Remez algorithm.
%   EXCHANGE(XK, H, METHOD, F, P, Q, W) performs one step of the Remez algorithm
%   for the best rational approximation of the CHEBFUN F of the target function
%   according to the first method (METHOD = 1), i.e. exchanges only one point,
%   or the second method (METHOD = 2), i.e. exchanges all the reference points.
%   XK is a column vector with the reference, H is the levelled error, P is the
%   numerator, and Q is the denominator of the trial
%   rational function P/Q and W is the weight function.
%
%   [XK, NORME, E_HANDLE, FLAG] = EXCHANGE(...) returns the modified reference
%   XK, the supremum norm of the error NORME (included as an output argument,
%   since it is readily computed in EXCHANGE and is used later in REMEZ), a
%   function handle E_HANDLE for the error, and a FLAG indicating whether there
%   were at least N+2 alternating extrema of the error to form the next
%   reference (FLAG = 1) or not (FLAG = 0).
%
%   [XK, ...] = EXCHANGE([], 0, METHOD, F, P, Q, N + 2) returns a grid of N + 2
%   points XK where the error F - P/Q alternates in sign (but not necessarily
%   equioscillates). This feature of EXCHANGE is useful to start REMEZ from an
%   initial trial function rather than an initial trial reference.

% Compute extrema of the error.
% Rational case:

rr = findExtrema(f, p, q, rh, h, xk);
disp('finished extrema search');
err_handle = @(x) feval(f, x) - rh(x);

% Select exchange method.
if ( method == 1 )                           % One-point exchange.
    [ignored, pos] = max(abs(feval(err_handle, rr)));
    pos = pos(1);
else                                           % Full exchange.
    pos = find(abs(err_handle(rr)) >= abs(h)); % Values above levelled error
end

% Add extrema nearest to those which are candidates for exchange to the
% existing exchange set.
[r, m] = sort([rr(pos) ; xk]);
v = ones(Npts, 1);
v(2:2:end) = -1;
er = [feval(err_handle, rr(pos)) ; v*h];
er = er(m);

% Delete repeated points.
repeated = diff(r) == 0;
r(repeated) = [];
er(repeated) = [];


% Determine points and values to be kept for the reference set.
s = r(1);    % Points to be kept.
es = er(1);  % Values to be kept.
for i = 2:length(r)
    if ( (sign(er(i)) == sign(es(end))) && (abs(er(i)) > abs(es(end))) )
        % Given adjacent points with the same sign, keep one with largest value.
        s(end) = r(i);
        es(end) = er(i);
    elseif ( sign(er(i)) ~= sign(es(end)) )
        % Keep points which alternate in sign.
        s = [s ; r(i)];    %#ok<AGROW>
        es = [es ; er(i)]; %#ok<AGROW>
    end
end

%xx = linspace(-1,1,100000);
%plot(xx, err_handle(xx));
%hold on
%plot(s, es,'or');
%hold off
%pause()

% Of the points we kept, choose n + 2 consecutive ones 
% that include the maximum of the error.
[norme, index] = max(abs(es));
d = max(index - Npts + 1, 1);
if ( Npts <= length(s) )
    %disp('Large candidate set:');
    %length(s)
    xk = s(d:d+Npts-1);
    flag = 1;
else
    xk = s;
    flag = 0;
end

end

function xk = makeNewReference(x0, x1)

xk = x0;
if ( length(x0) ~= length(x1) )
    error('makeNewReference:not equal points in both references')    
end

for i = 1:length(x0)
    % Find the index of the nearest point:
    [~, idx] = min(abs(x1-x0(i)));
    xk(i) = x0(i) + (x1(idx) - x0(i))/2;
end
end




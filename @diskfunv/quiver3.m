function varargout = quiver3( F, varargin )
%QUIVER3   3-D quiver plot of a DISKFUNV.
%   This file is here in case we want it, but we have not fully implemented
%   the use of three components in diskfunv. 
%
%   QUIVER3(F) plots velocity vectors as arrows with components F(1), F(2),
%   F(3), which are CHEBFUN2 objects. QUIVER3 automatically scales the arrows to
%   fit. The arrows are plotted on a uniform grid.
%
%   QUIVER3(Z,F) plots velocity vectors at the equally spaced surface points
%   specified by the matrix or CHEBFUN2 Z. If Z is a CHEBFUN2 then we use Z to
%   map the uniform grid.
%
%   QUIVER3(X,Y,Z,F) plots velocity vectors at (x,y,z). If X, Y, Z are CHEBFUN2
%   objects then we use X, Y, Z to map the uniform grid.
%
%   QUIVER3(F,S), QUIVER3(Z,F,S) or QUIVER3(X,Y,Z,F,S) automatically scales the
%   arrows to fit and then stretches them by S. Use S=0 to plot the arrows with
%   the automatic scaling.
%
%   QUIVER3(...,LINESPEC) uses the plot linestyle specified for the velocity
%   vectors.  Any marker in LINESPEC is drawn at the base instead of an arrow on
%   the tip.  Use a marker of '.' to specify no marker at all.  See PLOT for
%   other possibilities.
%
%   QUIVER(...,'numpts',N) plots arrows on a N by N uniform grid.
%
%   QUIVER3(...,'filled') fills any markers specified.
%
%   H = QUIVER3(...) returns a quiver object.
%
%   If F is a DISKFUNV with two components then we recommend using
%   DISKFUNV/QUIVER.
%
% See also QUIVER.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

numpts = 20; 

% Empty check:
if ( isempty( F ) )
    quiver([])
    return
end

if ( isempty(varargin) )
    varargin = {};
end

% Number of points to plot
j = 1; 
argin = {};
while ( ~isempty(varargin) )
    if strcmpi(varargin{1}, 'numpts')
        numpts = varargin{2};
        varargin(1:2) = [];
    else
        argin{j} = varargin{1};
        varargin(1) = [];
        j = j+1;
    end
end
varargin = argin; 

if ( isa(F, 'diskfunv') && ...
        (isempty(varargin) || ~isa(varargin{1}, 'diskfunv')) ) 
    % quiver3(F,...)
    
    nF = F.nComponents; 
    if ( nF == 2 ) 
        h = quiver( F, varargin{:} ); 
    else
        % Plot arrows are equally spaced locations
                [xx, yy] = diskpts(numpts, 0, 0);
                zz = zeros( size( xx )) ;
        vals1 = feval(F.components{1}, xx, yy, 'cart');
        vals2 = feval(F.components{2}, xx, yy, 'cart');
        vals3 = feval(F.components{3}, xx, yy, 'cart');
        h = quiver3( xx, yy, zz, vals1, vals2, vals3, varargin{:} );
        axis(1.1*[-1 1 -1 1])
    end

elseif ( isa(F, 'diskfun') && isa(varargin{1}, 'diskfunv') )  
    % quiver(Z,F,...)
    
    % Extract arguments:
    Z = F; 
    F = varargin{1};
    % Domain check:
    if ( ~domainCheck(Z, F.components{1} ) )
        error('CHEBFUN:DISKFUNV:quiver3:domain', ...
            'Object are not on the same domain.');
    end
    % Workout plotting locations:
    zz = new_data_locations( Z, numpts );
    %Plot: 
    h = quiver3( zz, F, varargin{2:end} );
    
elseif ( isa(F, 'double') && isa(varargin{1}, 'diskfunv') )  
    % quiver(zz,F,...)
    
    nF = F.nComponents; 
    if ( nF == 2 ) 
        h = quiver( F, varargin{:} ); 
    else
        % Plot arrows are equally spaced locations
        %dom = F.components{1}.domain;
        %x = linspace(dom(1), dom(2), numpts);
        %y = linspace(dom(3), dom(4), numpts);
        %[xx, yy] = meshgrid(x, y); 
        [xx,yy] = diskpts(numpts,0,0);
        vals1 = feval(F.components{1}, xx, yy, 'cart');
        vals2 = feval(F.components{2}, xx, yy, 'cart');
        vals3 = feval(F.components{3}, xx, yy, 'cart');
        h = quiver3( F, vals1, vals2, vals3, varargin{:} );
    end
    
elseif ( nargin > 3 )        
    % quiver(xx,yy,zz,F,...) or quiver(X,Y,Z,F,...)
    
    if ( isa(F,'double') )       
        % quiver(xx,yy,zz,F,...)
        
        % Extract arguments: 
        xx = F; 
        yy = varargin{1}; 
        zz = varargin{2}; 
        F = varargin{3};
        % Check that we have the right input types: 
        if ( ~isa(yy,'double') || ~isa(zz,'double') || ~isa(F,'diskfunv') )
            error('CHEBFUN:DISKFUNV:quiver3:badInputs1',...
                'Unrecognised input arguments.');
        end
        % Get plotting data: 
        dom = F.components{1}.domain;
        xdata = linspace( dom(1), dom(2), size(xx,1) );
        ydata = linspace (dom(3), dom(4), size(yy,2) ); 
        [mxdata, mydata]=meshgrid(xdata, ydata);
        [mxdata, mydata] = diskpts(0, mxdata, mydata);
        vals1 = feval(F.components{1}, mxdata, mydata, 'cart');
        vals2 = feval(F.components{2}, mxdata, mydata, 'cart');
        vals3 = feval(F.components{3}, mxdata, mydata, 'cart');
        % Plot:
        h = quiver3(xx, yy, zz, vals1, vals2, vals3, varargin{4:end});
        
    elseif isa(F,'diskfun')                         
        % quiver(X,Y,Z,F,...)
        
        % Extract arguments: 
        X = F; 
        Y = varargin{1}; 
        Z = varargin{2}; 
        F = varargin{3};
        % Check that we have the right input types: 
        if ( ~isa(Y,'diskfun') || ~isa(Z,'diskfun') || ~isa(F,'diskfunv') )
            error('CHEBFUN:DISKFUNV:quiver3:badInputs2', ...
                'Unrecognised input arguments.');
        end
        % Check domains: 
        if ( ~domainCheck(X, Y) || (~domainCheck(Y, Z) ) )
            error('CHEBFUN:DISKFUNV:quiver3:domain',...
                'Object are not on the same domain.');
        end
        
        % Get new data locations.
        xx = new_data_locations(X, numpts);
        yy = new_data_locations(Y, numpts);
        zz = new_data_locations(Z, numpts);
        
        % Plot quiver3.
        h = quiver3( xx, yy, zz, F, varargin{4:end} );
    end
    
else
    error('CHEBFUN:DISKFUNV:quiver3:badInputs3',...
        'Unrecognised input arguments.');
end

if ( nargout > 0 )
    varargout = { h };
end



function newloc = new_data_locations( f1 , numpts )
% Generate new arrow location if first two inputs are CHEBFUN2 objects.

% Check the CHEBFUN2 objects are on the same domain.
%dom = f1.domain;
%dom = [-1 1 -1 1];
% mesh 'em up for the quiver arrows.
x = linspace(-1, 1, numpts);
y = linspace(-1, 1, numpts);

[a, b] = meshgrid(x, y);
[xx, yy] = diskpts(0, a, b); 
% Use CHEBFUN2 to generate data locations.
newloc = feval(f1, xx, yy,'cart');      
end

function [xx,yy] = diskpts(s, x,y)  % gets good spaced points on disk in cartesian coords. Will need improvement.
                                    %for user-input coords, this throws out
                                    % choices that are not on the unit
                                    % disk. 
if (s > 0) 
x = linspace(-1,1,s); 
y = linspace(-1,1,s);
[xx, yy] = meshgrid(x,y); 
else
xx = x; 
yy = y;
end

xx=xx(:);
yy = yy(:);

R = xx.^2+yy.^2;   
a = find(bsxfun(@max, R, ones(length(xx),1))-R); %throw out stuff outside unit disk                                           
xx = xx(a);
yy = yy(a);

end

end





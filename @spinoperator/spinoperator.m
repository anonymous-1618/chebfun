classdef spinoperator
%SPINOPERATOR   Abstract class for representing the spatial part of differential 
%operators for time-dependent PDEs.
%
% See also SPINOP, SPINOP2, SPINOP3.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = public )
        dimension           % Spatial dimension
        domain              % Spatial domain
        linearPart          % Linear part of the operator
        nonlinearPart       % Nonlinear part of the operator
        numVars             % Number of unknown functions (>1 for systems)
    end
    
    properties ( Access = public, Hidden = true )
        nonlinearPartCoeffs % Nonlinear part of the operator in coeffs space
        nonlinearPartVals   % Nonlinear part of the operator in values space
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% ABSTRACT AND NON-STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Abstract = true, Static = false )
        
        % Discretize a SPINOPERATOR:
        [L, Nc] = discretize(S, N)
        
        % Initialize a movie when solving a PDE specified by a SPINOPERATOR:
        [p, plotOption] = initializeMovie(S, dt, pref, v, gridPoints)
        
        % Plot a movie when solving a PDE specified by a SPINOPERATOR:
        plotMovie
        
    end
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CONCRETE AND STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Abstract = false, Static = true )
        
        [uout, tout] = solvepde(varargin)
   
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CONCRETE AND NON-STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Abstract = false, Static = false )
        
        % Create a contour around each eigenvalue of the linear part of a 
        % SPINOPERATOR:
        LR = computeLR(S, dt, L, M, N)
        
        % Get the nonlinear parts, in coefficient and value space, of a
        % SPINOPERATOR:
        [Nc, Nv] = getNonlinearPartsCoeffsAndVals(S)
           
    end
    
end
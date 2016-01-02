classdef spinpref2 < spinpreference
%SPINPRE2   Class for managing SPIN2/SPINOP2 preferences.
%
% Available preferences ([] = defaults):
%
%   dealias                   * If 1, use the 2/3-rule to zero high wavenumbers.
%     [0]                       No dealiasing by default.
% 
%   dt                        * Time-step for time discretization. Default is 
%     []                        empty, i.e., adaptively chosen by the code to 
%                               achieve errTol. 
%
%   dtmax                     * Maximum time-step when using an apative grid in
%     [1]                       time.
%
%   dtmin                     * Minimum time-step when using an apative grid in
%     [1e-10]                   time.
%
%   errTol                    * Desired accuracy on the solution.
%     [1e-4]
%
%   M                         * Number of points for complex means to evaluate
%     [64]                      the phi-functions.
%
%   N                         * Number points in each direction for spatial 
%     []                        discretization. Default is empty, i.e., 
%                               adaptively chosen by the code to achieve errTol.
%
%   Nmin                      * Minimum number of points in each direction when 
%     [128]                     using an adaptive grid in space.
%
%   Nmax                      * Maximum number of points in each direction when   
%     [1024]                    using an adaptive grid in space.
%                                         
%   plotting                  * Plotting options: 'movie' for plotting a 
%     ['movie']                 movie of the solution, [] for no plotting.
%
%   scheme                    * Time-stepping scheme.
%     [@etdrk4]  
%      @exprk5s8
%      @krogstad
%      @eglm433
%      @pecec433
%
%   view                      * Graph view specification when using 'movie'.
%     [0 90]   
%
% See also SPINPREF, SPINPREF3.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = public )
        view = [0 90];        % Graph view specification
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function pref = spinpref2(varargin) 
            if ( nargin == 0 )
                pref.dtmax = 1;
                pref.errTol = 1e-4;
                pref.iterPlot = 4;
                pref.Nmin = 128;
                pref.Nmax = 1024;
            else
                pref = spinpref2();
                for k = 1:nargin/2
                    pref.(varargin{2*(k-1)+1}) = varargin{2*k};
                end
            end
        end
    end
    
end
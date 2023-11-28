function OPTION=param(OPTION)
%
% Create SDPA parameters. 
% If there is no argument, default parameters are returned.
% 					
% OPTION=param         % for default parameter
% or 
% OPTION=param(field1,value1,field2,value2,....)
% 
% <INPUT>
% - field?: string : field name
% - value?: numeric or string : 
%
% <OUTPUT>
% - OPTION: structure data: each field is as follows:
% * maxIteration : The maximum number of iterations. 
% * epsilonStar  : The accuracy of an approximate optimal solution
%                  for primal and dual SDP.
% * lambdaStar   : An initial point.
% * omegaStar    : The search region for an optimal solution.
% * lowerBound   : Lower bound of the minimum objective value of 
%                  the primal SDP.
% * upperBound   : Upper bound of the maximum objective value of 
%                  the dual SDP 
% * betaStar     : The parameter for controlling the search direction
%                  if the current point is feasible.
% * betaBar      : The parameter for controlling the search direction
%                  if the current point is infeasible.
% * gammaStar    : A reduction factor for the primal and dual step 
%                  lengths.
% * epsilonDash  : The relative accuracy of an approximate optimal
%                  solution between primal and dual SDP.
% * isSymmetric  : The flag for the checking the symmetricity of input 
%                  matrices. (0 => no check, 1=> check)
% * isDimacs     : The flag to compute DIMACS ERROR
%                  (0 => no computation, 1=> computation)
% * xPrint        : (default %+8.3e, NOPRINT skips printout)
% * XPrint        : (default %+8.3e, NOPRINT skips printout)
% * YPrint        : (default %+8.3e, NOPRINT skips printout)
% * infPrint      : (default %+10.16e, NOPRINT skips printout)
% * print         : Destination of file output. the default setting is 
%                   stdout by 'display'.
%                   If print is set 'no' or empty, no message 
%                   is print out
% * resultFile    : Destination of detail file output
% * NumThreads    : Number of Threads for internal computation

% This file is a component of SDPA
% Copyright (C) 2004-2012 SDPA Project
% 
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
% 
% SDPA-M: $Revision: 6.2 $
% $Id: param.m,v 6.2 2005/05/28 02:36:40 drophead Exp $

% create default OPTION
OPTION0.maxIteration = 100;
OPTION0.epsilonStar  = 1.0E-7;
OPTION0.lambdaStar   = 1.0E2;
OPTION0.omegaStar    = 2.0;
OPTION0.lowerBound   = -1.0E5;
OPTION0.upperBound   = 1.0E5;
OPTION0.betaStar     = 0.1;
OPTION0.betaBar      = 0.2;
OPTION0.gammaStar    = 0.9;
OPTION0.epsilonDash  = 1.0E-7;
OPTION0.isSymmetric  = 0;
OPTION0.isDimacs     = 0;
OPTION0.xPrint       = '%+8.3e';
OPTION0.XPrint       = '%+8.3e';
OPTION0.YPrint       = '%+8.3e';
OPTION0.infPrint     = '%+16.10e';
OPTION0.print        = 'display';
OPTION0.resultFile   = '';
try 
    OPTION0.NumThreads   = maxNumCompThreads; % Max Avialable Number
catch
    fprintf(['Function maxNumCompThreads is not found, NumThreads ' ...
             'is set as 1.\n']);
    OPTION0.NumThreads   = 1;
end

    
OPTION0.aggConeSize  = [];

if (nargin == 0) || isempty(OPTION)
    OPTION = OPTION0;
    return
else
    if ~isfield(OPTION,'maxIteration')
        OPTION.maxIteration=OPTION0.maxIteration;
    elseif ~isnumeric(OPTION.maxIteration)
        error('OPTION.maxIteration must be numeric.');
    end
    %
    if ~isfield(OPTION,'epsilonStar')
        OPTION.epsilonStar=OPTION0.epsilonStar;
    elseif ~isnumeric(OPTION.epsilonStar)
        error('epsilonStar must be numeric.');
    end
    %
    if ~isfield(OPTION,'lambdaStar')
        OPTION.lambdaStar=OPTION0.lambdaStar;
    elseif ~isnumeric(OPTION.lambdaStar)
        error('OPTION.lambdaStar must be numeric.');
    end
    %
    if ~isfield(OPTION,'omegaStar')
        OPTION.omegaStar=OPTION0.omegaStar;
    elseif ~isnumeric(OPTION.omegaStar)
        error('OPTION.omegaStar must be numeric.');
    end
    %
    if ~isfield(OPTION,'lowerBound')
        OPTION.lowerBound=OPTION0.lowerBound;
    elseif ~isnumeric(OPTION.lowerBound)
        error('OPTION.lowerBound must be numeric.');
    end
    %
    if ~isfield(OPTION,'upperBound')
        OPTION.upperBound=OPTION0.upperBound;
    elseif ~isnumeric(OPTION.upperBound)
        error('OPTION.upperBound must be numeric.');
    end
    %
    if ~isfield(OPTION,'betaStar')
        OPTION.betaStar=OPTION0.betaStar;
    elseif ~isnumeric(OPTION.betaStar)
        error('OPTION.beaStar must be numeric.');
    end
    %
    if ~isfield(OPTION,'betaBar')
        OPTION.betaBar=OPTION0.betaBar;
    elseif ~isnumeric(OPTION.betaBar)
        error('OPTION.betaBar must be numeric.');
    end
    %
    if ~isfield(OPTION,'gammaStar')
        OPTION.gammaStar=OPTION0.gammaStar;
    elseif ~isnumeric(OPTION.gammaStar)
        error('OPTION.gammaStar must be numeric.');
    end
    %
    if ~isfield(OPTION,'epsilonDash')
        OPTION.epsilonDash=OPTION0.epsilonDash;
    elseif ~isnumeric(OPTION.epsilonDash)
        error('OPTION.epsilonDash must be numeric.');
    end
    %
    if isfield(OPTION,'searchDir')
        disp('Parameter *searchDir* is no longer supported.');
        disp('HRVW/KSH/M is automatically used.');
    end
    %
    if ~isfield(OPTION,'isSymmetric')
        OPTION.isSymmetric=OPTION0.isSymmetric;
    elseif ~isnumeric(OPTION.isSymmetric) || ((OPTION.isSymmetric~=0) && (OPTION.isSymmetric~=1))
        error('OPTION.isSymmetric must be 0 or 1.');
    end
    %
    if ~isfield(OPTION,'isDimacs')
        OPTION.isDimacs=OPTION0.isDimacs;
    elseif ~isnumeric(OPTION.isDimacs) || ((OPTION.isDimacs~=0) && (OPTION.isDimacs~=1))
        error('OPTION.isDimacs must be 0 or 1.');
    end
    %
    if ~isfield(OPTION,'XPrint')
        OPTION.XPrint=OPTION0.XPrint;
    elseif ~ischar(OPTION.XPrint)
        error('OPTION.XPrint must be string.');
    end
    %
    if ~isfield(OPTION,'YPrint')
        OPTION.YPrint=OPTION0.YPrint;
    elseif ~ischar(OPTION.YPrint)
        error('OPTION.YPrint must be string for printf.');
    end
    %
    if ~isfield(OPTION,'infPrint')
        OPTION.infPrint=OPTION0.infPrint;
    elseif ~ischar(OPTION.infPrint)
        error('OPTION.infPrint must be string for printf.');
    end
    %
    if isfield(OPTION,'print') && ...
            (isempty(OPTION.print) || length(OPTION.print) == 0)
        OPTION.print = 'no';
    end
    if ~isfield(OPTION,'print')
        OPTION.print=OPTION0.print;
    elseif ~ischar(OPTION.print)
        disp('*** OPTION.print must be string for FILE. ***');
        disp('  "display" is for stdout.');
        disp('  "no" or empty is for skip message.');
        disp('  filename is filename in which message will be written.');
        error('*** OPTION.print must be string for FILE. ***');
    end
    %
    if ~isfield(OPTION,'resultFile') || isempty(OPTION.resultFile)
        OPTION.resultFile=OPTION0.resultFile;
    elseif ~ischar(OPTION.resultFile)
        error('OPTION.resultFile must be string.');
    end
    %
    if ~isfield(OPTION,'NumThreads')
        OPTION.NumThreads=OPTION0.NumThreads;
    elseif ~isnumeric(OPTION.NumThreads)
        error('OPTION.NumThreads must be positive integer.');
    end
    
    if ~isfield(OPTION,'aggConeSize')
        OPTION.aggConeSize = OPTION0.aggConeSize;
    elseif (~isempty(OPTION.aggConeSize)) && ...
            ((~isnumeric(OPTION.aggConeSize)) || (OPTION.aggConeSize <=0))
        error('OPTION.aggConeSize must be a positive integer.');
    end
end

return

% End of File

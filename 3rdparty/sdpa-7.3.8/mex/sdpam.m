function [objVal,x,X,Y,INFO]=sdpam(mDIM,nBLOCK,bLOCKsTRUCT,c,F,...
                                   x0,X0,Y0,OPTION)
%
% Compute the solution of standard SDP.
% Since some of input arguments are optional, sdpam can be
% overloaded as below.
%               
% [objVal,x,X,Y,INFO] = sdpam(mDIM,nBLOCK,bLOCKsTRUCT,c,F,
%                             x0,X0,Y0,OPTION);
%
% <INPUT>
% - mDIM       : integer   ; number of primal variables
% - nBLOCK     : integer   ; number of blocks of F
% - bLOCKsTRUCT: vector    ; represetns the block structure of F
% - c          : vector    ; coefficient vector
% - F          : cell array; coefficient matrices
% - x0,X0,Y0   : cell array; initial point
% - OPTION     : structure ; options
% 
% <OUTPUT>
% - objVal: [objValP objValD]; optimal value of P and D
% - x     : vector           ; optimal solution
% - X,Y   : cell arrray      ; optimal solutions
% - INFO  : structure        ; infomation of the solution
% 

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
% $Id: sdpam.m,v 6.2 2005/05/28 02:36:40 drophead Exp $

t = cputime; 

if (nargin < 5 | nargin > 9)
  error('incorrect number of input arguments')
elseif nargin == 5
  % make initial points empty
  x0=[];X0=[];Y0=[];
  % load default parameters
  OPTION=param;
  % solve by SDPA
  [objVal,x,X,Y,INFO]=mexsdpa(mDIM,nBLOCK,bLOCKsTRUCT,...
                                 c,F,x0,X0,Y0,OPTION);
elseif nargin == 6
  % use OPTION given by arguments
  OPTION=param(x0);
  % make initial points empty
  x0=[];X0=[];Y0=[];
  [objVal,x,X,Y,INFO]=mexsdpa(mDIM,nBLOCK,bLOCKsTRUCT,...
                                 c,F,x0,X0,Y0,OPTION);
elseif nargin == 8
  % load default parameters
  OPTION=param;
  %solve by SDPA
  [objVal,x,X,Y,INFO]=mexsdpa(mDIM,nBLOCK,bLOCKsTRUCT,...
                                 c,F,x0,X0,Y0,OPTION);
elseif nargin == 9
  OPTION=param(OPTION);
  % solve by SDPA
  [objVal,x,X,Y,INFO]=mexsdpa(mDIM,nBLOCK,bLOCKsTRUCT,...
                                 c,F,x0,X0,Y0,OPTION);
end
INFO.cpusec = cputime-t;

% End of File

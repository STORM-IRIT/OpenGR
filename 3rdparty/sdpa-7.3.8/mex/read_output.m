function [objVal,x,X,Y,INFO] = read_output(filename,mDIM,nBLOCK,bLOCKsTRUCT);
%
% Read results from SDPA output file
%
% [objVal,x,X,Y,INFO] = read_output(filename,mDIM,nBLOCK,bLOCKsTRUCT);
%
% <INPUT>
% - filename: string; filename of the SDPA output
% - mDIM       : integer; number of primal variables
% - nBLOCK     : integer; number of blocks of F
% - bLOCKsTRUCT: vector; represetns the block structure of F
%
% <OUTPUT>
% - objVal: [objValP objValD]; optimal value of P and D
% - x     : vector           ; optimal solution
% - X,Y   : cell arrray      ; optimal solutions
% - INFO  : structure        ; infomation of the solution
%    
% /* -------------------------------------------------------------
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
% ------------------------------------------------------------- */
     
if (nargin~=4)
    error('incorrect number of input arguments');
end

if max(size(bLOCKsTRUCT)) ~= nBLOCK
    error('Inconsinstent between nBLOCK and bLOCKsTRUCT');
end
% for filename starting with '~'
if filename(1) == '~'
    filename = strcat(getenv('HOME'),filename(2:length(filename)));
end

[objVal,x,X,Y,INFO] = ...
    mexReadOutput(filename,full(mDIM),full(nBLOCK),full(bLOCKsTRUCT));
    
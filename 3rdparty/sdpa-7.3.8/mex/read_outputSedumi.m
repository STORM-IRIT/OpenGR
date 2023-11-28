function [x,y,info] = read_outputSedumi(filename,m,K);
% [x,y,info] = read_outputSedumi(filename,m,K);
%
% <INPUT>
% - filename: string; filename of the SDPA output
% - m       : integer; number of dual variables (that is, length(b))
% - K       : struct; cone information of SeDuMi
%
% <OUTPUT>
% - x,y   : vector           ; optimal solutions
% - info  : structure        ; infomation of the solution
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
 
nBlock = 0;
blockStruct = [];
if (~isfield(K,'l')) || (isempty(K.l))
    K.l = 0;
else
    nBlock = 1;
    blockStruct = [-K.l];
end

if isfield(K,'s') && ~isempty(K.s)
    if size(K.s,2) ~= 1
        % fprintf('Transposing K.s to a column vector\n');
        K.s = K.s';
    end
    nBlock = nBlock + length(K.s);
    blockStruct = [blockStruct; K.s];
end

[objValo,xo,Xo,Yo,INFOo] = read_output(filename,m,nBlock, ...
                                       blockStruct);

if isempty(INFOo) == 1
    fprintf('Output [x,y,info] will be empty');
    x = [];
    y = [];
    info = [];
end

    
info.primal = -objValo(1);
info.dual   = -objValo(2);
info.phaseValue = INFOo.phasevalue;
info.iter   = INFOo.iteration;
info.cpusec = INFOo.cpusec;

y = -xo;

x = [];
coneIndex = 1;
if K.l > 0
    x = [x; reshape(Yo{1},K.l,1)];
    coneIndex = coneIndex + 1;
end
while coneIndex <= length(blockStruct)
    blk = blockStruct(coneIndex);
    x = [x; reshape(Yo{coneIndex},blk*blk,1)];
    coneIndex = coneIndex + 1;
end

    

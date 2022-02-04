function [At,b,c,K,blockStruct] = SDPAToSedumi(filename);
% [At,b,c,K,blockStruct] = SDPAToSedumi(filename);
% [At,b,c,K,blockStruct] = SDPAToSedumi('example1.dat-s');
% A converter from SDPA sparse format to SeDuMi Input
% (C) SDPA Project 2009
%
% filename : Filename for SDPA dat-s
% A,b,c,K  : SeDuMi Input
% blockStruct : SDPA blockStruct information
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
    
    % for filename starting with '~'
    if filename(1) == '~'
        filename = strcat(getenv('HOME'),filename(2:length(filename)));
    end
    
    DEBUG_SWITCH = 0;
    
    if DEBUG_SWITCH == 1
        tic
    end
    [At,b,c,K,blockStruct] = mexReadSDPAToSedumi(filename);
    if DEBUG_SWITCH
        toc
    end

    if DEBUG_SWITCH == 1
       tic
           [At2,b2,c2,K2] = fromsdpa(filename);
       toc
           fprintf('diff A[%e], b[%e], c[%e], K.l[%d], K.s[%d]\n',...
                   nnz(At-At2), norm(b-b2), norm(c-c2), K.l-K2.l,...
                   sum(abs(K.s'-K2.s)));
    end



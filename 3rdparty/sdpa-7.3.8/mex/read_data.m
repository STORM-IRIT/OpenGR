function [mDIM,nBLOCK,bLOCKsTRUCT,c,F]=read_data(filename); 
%
% Read a problem in SDPA sparse format.
%
% [mDIM,nBLOCK,bLOCKsTRUCT,c,F] = read_data(fname)
%
% <INPUT>
% - filename: string; filename of the SDP data with SDPA foramt.
%
% <OUTPUT>
% - mDIM       : integer; number of primal variables
% - nBLOCK     : integer; number of blocks of F
% - bLOCKsTRUCT: vector; represetns the block structure of F
% - c          : vector; coefficient vector
% - F          : cell array; coefficient matrices
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
% $Id: read_data.m,v 6.2 2005/05/28 02:36:40 drophead Exp $

% check the validity of the arguments
if ( nargin ~= 1 | ( nargin == 1 & ~ischar(filename) ) )
  error('input argument must be a filename');
end

% identify whether a file is sparse format or not.
bsparse=0;
len=length(filename);
if len >= 2
  str=filename(end-1:end);
  if strncmp(str,'-s',2) 
    bsparse=1;
  end
end

fid=fopen(filename,'r');
if fid == -1
  error(sprintf('Cannot open %s',filename));
end

% skip comment and after it, read a number of decision variables (mDIM)
while 1
  str=fgetl(fid);
  if( str(1)~='*' & str(1) ~='"' )
    mDIM=sscanf(str,'%d',1);
    break;
  end
end
%disp(sprintf('mDIM=%d',mDIM));

% read a number of blocks (nBLOCK)
nBLOCK=fscanf(fid,'%d',1);
%disp(sprintf('nBLOCK=%d',nBLOCK));

% read each size of blocks (bLOCKsTRUCT)
bLOCKsTRUCT=zeros(nBLOCK,1);
for idx=1:nBLOCK
  bLOCKsTRUCT(idx)=fscanf(fid,'%*[^0-9+-]%d',1);
  if bLOCKsTRUCT(idx) == 1
    bLOCKsTRUCT(idx) = -1;
  end
end

% read cost vector (c)
c=zeros(mDIM,1);
for idx=1:mDIM
  c(idx)=fscanf(fid,'%*[^0-9+-]%lg',1);
end

% read coefficient matrices (F)
F=cell(nBLOCK,mDIM+1);

if bsparse
  % sparse format case
  while 1
    [k,cnt]=fscanf(fid,'%*[^0-9+-]%d',1);
    [l,cnt]=fscanf(fid,'%*[^0-9+-]%d',1);
    [i,cnt]=fscanf(fid,'%*[^0-9+-]%d',1);
    [j,cnt]=fscanf(fid,'%*[^0-9+-]%d',1);
    [value,cnt]=fscanf(fid,'%*[^0-9+-]%lg',1);
    if cnt
      if isempty(F{l,k+1})
	size=abs(bLOCKsTRUCT(l));
	if bLOCKsTRUCT(l) < 0
	  F{l,k+1}=sparse(zeros(size,1));
	else
	  F{l,k+1}=sparse(zeros(size));
	end
      end
      if bLOCKsTRUCT(l) < 0
	F{l,k+1}(i)=value;
      else
        if j < i
          tmpj = j;
          j = i;
          i = tmpj;
        end 
	if i < j
	  F{l,k+1}(i,j)=value;
	  F{l,k+1}(j,i)=value;
	elseif i == j
	  F{l,k+1}(i,j)=value;
	end
      end
    else 
      break;
    end
  end
else
  % dense format case
  for k=1:mDIM+1
    for l=1:nBLOCK
      size=abs(bLOCKsTRUCT(l));
      if bLOCKsTRUCT(l) > 0
	F{l,k}=zeros(size);
	for i=1:size
	  for j=1:size
	    [value,cnt]=fscanf(fid,'%*[^0-9+-]%lg',1);
	    if cnt
	      F{l,k}(i,j)=value;
	    else
	      error(sprintf('Failed to read an element at %d %d %d %d'),...
		k-1,l,i,j);
	    end
	  end
	end
      else
	F{l,k}=zeros(size,1);
	for i=1:size
	  [value,cnt]=fscanf(fid,'%*[^0-9+-]%lg',1);
	  if cnt
	    F{l,k}(i)=value;
	  else
	    error(sprintf('Failed to read an element at %d %d %d %d'),...
	      k-1,l,i,i);
	  end
	end
      end
    end
  end
end
fclose(fid);

% End of File


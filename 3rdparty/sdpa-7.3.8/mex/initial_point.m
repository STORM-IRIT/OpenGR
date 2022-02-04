function [x0,X0,Y0] = initial_point(filename,mDIM,nBLOCK,bLOCKsTRUCT)
%
% Read in a sparse initial point file
%
% [x0,X0,Y0] = initial_point(fname, mDIM, nBLOCK, bLOCKsTRUCT)
%
% <INPUT>
% - filename   : string ; filename of initial point file
% - mDIM       : integer; number of primal variables
% - nBLOCK     : integer; number of blocks of F
% - bLOCKsTRUCT: vector ; represetns the block structure of F
%
% <OUTPUT>
% - x0: vector    ; initial point of x
% - X0: cell array; initial point of X
% - Y0: cell array; initial point of Y
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
% $Id: initial_point.m,v 6.2 2005/05/28 02:36:40 drophead Exp $

% check the validity of arguments
if nargin ~= 4
  error('input arguments must be 4.');
else
  if ~ischar(filename)
    error('1st argument must be a filename.');
  end
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
% disp(sprintf('sparse=%d',bsparse));

fid=fopen(filename,'r');
if fid == -1
  error(sprintf('Cannot open %s',filename));
end

x0=zeros(mDIM,1);
X0=cell(nBLOCK,1);
Y0=cell(nBLOCK,1);

% read initial point x0
for idx=1:mDIM
  x0(idx)=fscanf(fid,'%*[^0-9+-]%lg',1);
end

% read initial points X0, Y0
if bsparse
  % sparse format case
  while 1
    [k,cnt]=fscanf(fid,'%*[^0-9+-]%d',1);
    [l,cnt]=fscanf(fid,'%*[^0-9+-]%d',1);
    [i,cnt]=fscanf(fid,'%*[^0-9+-]%d',1);
    [j,cnt]=fscanf(fid,'%*[^0-9+-]%d',1);
    [value,cnt]=fscanf(fid,'%*[^0-9+-]%lg',1);
    if cnt
      if k == 1
	% X0
	if isempty(X0{l})
	  size=abs(bLOCKsTRUCT(l));
	  if bLOCKsTRUCT(l) < 0
	    X0{l}=sparse(zeros(size,1));
	  else
	    X0{l}=sparse(zeros(size));
	  end
	end
	if bLOCKsTRUCT(l) < 0
	  X0{l}(i)=value;
	else
	  if i < j
	    X0{l}(i,j)=value;
	    X0{l}(j,i)=value;
	  elseif i == j
	    X0{l}(i,j)=value;
	  end
	end
      elseif k==2
	% Y0
	if isempty(Y0{l})
	  size=abs(bLOCKsTRUCT(l));
	  if bLOCKsTRUCT(l) < 0
	    Y0{l}=sparse(zeros(size,1));
	  else
	    Y0{l}=sparse(zeros(size));
	  end
	end
	if bLOCKsTRUCT(l) < 0
	  Y0{l}(i)=value;
	else
	  if i < j
	    Y0{l}(i,j)=value;
	    Y0{l}(j,i)=value;
	  elseif i == j
	    Y0{l}(i,j)=value;
	  end
	end
      end
    else 
      break;
    end
  end
else
  % dense format case
  
  % X0
  for l=1:nBLOCK
    size=abs(bLOCKsTRUCT(l));
    if bLOCKsTRUCT(l) > 0
      X0{l}=zeros(size);
      for i=1:size
	for j=1:size
	  [value,cnt]=fscanf(fid,'%*[^0-9+-]%lg',1);
	  if cnt
	    X0{l}(i,j)=value;
	  else
	    error(sprintf('Failed to read an element X0 at %d %d %d'),...
	      l,i,j);
	  end
	end
      end
    else
      X0{l}=zeros(size,1);
      for i=1:size
	[value,cnt]=fscanf(fid,'%*[^0-9+-]%lg',1);
	if cnt
	  X0{l}(i)=value;
	else
	  error(sprintf('Failed to read an element X0 at %d %d %d'),...
	    l,i,i);
	end
      end
    end
  end
  
  % Y0
  for l=1:nBLOCK
    size=abs(bLOCKsTRUCT(l));
    if bLOCKsTRUCT(l) > 0
      Y0{l}=zeros(size);
      for i=1:size
	for j=1:size
	  [value,cnt]=fscanf(fid,'%*[^0-9+-]%lg',1);
	  if cnt
	    Y0{l}(i,j)=value;
	  else
	    error(sprintf('Failed to read an element Y0 at %d %d %d'),...
	      l,i,j);
	  end
	end
      end
    else
      Y0{l}=zeros(size,1);
      for i=1:size
	[value,cnt]=fscanf(fid,'%*[^0-9+-]%lg',1);
	if cnt
	  Y0{l}(i)=value;
	else
	  error(sprintf('Failed to read an element Y0 at %d %d %d'),...
	    l,i,i);
	end
      end
    end
  end
end
fclose(fid);

% End of File

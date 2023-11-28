function [x,y,info]=sedumiwrap(A,b,c,K,pars,OPTION);

%
% SeDuMi wrapper for SDPA
%               
% [x,y,info]=sedumiwrap(A,b,c,K,pars,OPTION);
% or
% [x,y,info]=sedumiwrap(A,b,c,K);   % with SDPA-M default parameter
%
% Note :
% 'A', in each SDP block, only upper triangle part is used.
% 'K' can include only 'f'(free) 'l'(linear) 's'(SDP) cones.
% 'pars' information is NOT used (just for SeDuMi compatibility)
% 'OPTION' is option structure for SDPA-M (for details, try 'help param')
% 'info' information is based on SDPA-M
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
% SDPA-M: $Revision: 7.3 $

t = cputime; 

fprintf('-SeDuMi Wrapper for SDPA Start-\n');

if (nargin < 4 | nargin > 6)
    error('incorrect number of input arguments')
end
if nargin >= 5
    if isfield(OPTION,'print') && ~isempty(OPTION.print)
        fprintf('Note: pars information [5th argument] is not used\n');
    end
end
if nargin < 6
    OPTION = param;
else
    OPTION = param(OPTION);
end


if isfield(K,'q') && ~isempty(K.q)
   error('Current Wrapper cannot handle K.q'); 
end

if isfield(K,'r') && ~isempty(K.r)
   error('Current Wrapper cannot handle K.r'); 
end

if size(b,2) ~= 1
    % fprintf('Transposing b to a column vector');
    b = b';
end
if size(b,2) ~= 1
    error('b must be a vector');
end

if size(c,2) ~= 1
    % fprintf('Transposing c to a column vector\n');
    c = c';
end
if size(c,2) ~= 1
    error('c must be a vector');
end

%%%%%
% Constructing at least one SDP cone if necesary 
if (~isfield(K,'s')) || (isempty(K.s))
    [K] = LPtoLP_SDP(K); 
end
%%%%%

totalLength = 0;
Kf = 0;
if isfield(K,'f') && ~isempty(K.f)
    totalLength = totalLength + K.f;
    Kf = K.f;
end
if isfield(K,'l') && ~isempty(K.l)
    totalLength = totalLength + K.l;
else
    K.l = 0;
end
if isfield(K,'s') && ~isempty(K.s)
    if size(K.s,2) ~= 1
        % fprintf('Transposing K.s to a column vector\n');
        K.s = K.s';
    end
    Ks = sum(K.s .* K.s);
    totalLength = totalLength + Ks;
else
    error('Cannot handle empty K.s');
end
m = size(b,1);
n = size(c,1);

[mA,nA] = size(A);
if (m~=mA || n~=nA) && (m~=nA || n~=mA)
    fprintf('Inconsistent Size of A,b,c\n');
    fprintf('size(A) = [%d,%d], size(b) = %d, size(c) = %d\n',...
            mA,nA, m,n);
    error('Cannot continue...');
end

if (n~=totalLength)
    fprintf('Inconsistent Size of c and K\n');
    fprintf('size(c) = %d, totalSize(K) = %d\n',...
            n, totalLength);
    error('Cannot continue...');
end

if ~issparse(A)
    if isfield(OPTION,'print') && ~isempty(OPTION.print)
        fprintf('Converting A from dense to sparse\n');
    end
    A = sparse(A);
end

if issparse(b)
    % fprintf('Converting b from sparse to dense\n');
    b = full(b);
end

if ~issparse(c)
    % fprintf('Converting c from dense to sparse\n');
    c = sparse(c);
end

if issparse(K.s)
    % fprintf('Converting K.s from sparse to dense\n');
    K.s = full(K.s);
end

if Kf ~= 0
    if isfield(OPTION,'print') && ~isempty(OPTION.print)
        fprintf(['Free Variables are divided into positive and ' ...
                 'negative part of LP cone\n']);
    end
    Af = A(:,1:Kf);
    Kl = K.l;
    Al = A(:,Kf+1:Kf+Kl);
    As = A(:,Kf+Kl+1:Kf+Kl+Ks);
    Anew = [Af, -Af, Al, As];
        
    cf = c(1:Kf);
    cl = c(Kf+1:Kf+Kl);
    cs = c(Kf+Kl+1:Kf+Kl+Ks);
    cnew = [cf; -cf; cl; cs];
    
    Knew.l = 2*Kf + Kl;
    Knew.s = K.s;
    
    A = Anew;
    c = cnew;
    K = Knew;
    
    clear Af;
    clear Al;
    clear As;
    clear Anew;
    clear cf;
    clear cl;
    clear cs;
    clear cnew;
    clear Knew;
end

if isfield(K,'s') && ~isempty(K.s)
    if size(K.s,2) ~= 1
        K.s = Ks';
    end
end

%%%%%
% Aggregating small SDP cones into larger SDP cones
aggSW = 0;         
minNoSDPcones = 3; 
if (isfield(OPTION,'aggConeSize'))  && (~isempty(OPTION.aggConeSize)) && (isnumeric(OPTION.aggConeSize)) && ... 
        ((OPTION.aggConeSize > 0)) && (length(K.s') > minNoSDPcones)  && (length(find(K.s' < OPTION.aggConeSize)) > minNoSDPcones) 
    fprintf('OPTION.aggConeSize = %d\n',OPTION.aggConeSize)
    K0 = K; 
    aggSW = 1;  
    [A,c,K1] = aggSDPcones(A,c,K,OPTION.aggConeSize); 
    A = sparse(A);
    c = sparse(c); 
    K = K1; 
end
%%%%%

  
% A should be transposed when passed to mex
[mA,nA] = size(A);
if mA ~= K.l + sum(K.s.*K.s)
    A = A';
end
% fprintf('size(A) = (%d,%d)\n',size(A,2),size(A,1));
% fprintf('length(K.s) = %d\n',length(K.s)); 

[x,y,info] = mexSedumiWrap(A,b,c,K,OPTION);

%%%%%
% Retrieving the origianl primal SDP cone variables 
if aggSW == 1
    if (isfield(K,'l')) && (~isempty(K.l)) && (K.l > 0)
        xSDP = x(K.l+1:size(x,1),1);
        xLP = x(1:K.l,1);
        [xSDP] = mexDisAggSDPsol(xSDP,K0.s,K1.s); 
        x = [xLP; xSDP]; 
    else
        [x] = mexDisAggSDPsol(x,K0.s,K1.s); 
    end
end
%%%%%

if Kf ~=0
    xlength = size(x);
    xnew = x(1:Kf) - x(Kf+1:Kf+Kf);
    xnew = [xnew; x(Kf+Kf+1:xlength)];
    x = xnew;
end

info.cpusec = cputime-t;
% if isfield(OPTION,'print') && ~isempty(OPTION.print)
    fprintf('-SeDuMi Wrapper for SDPA End-\n');
% end

% End of File

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [K] = LPtoLP_SDP(K)

if isfield(K,'q') && ~isempty(K.q)
   error('Current Wrapper cannot handle K.q'); 
end

if isfield(K,'r') && ~isempty(K.r)
   error('Current Wrapper cannot handle K.r'); 
end

if isfield(K,'s') && ~isempty(K.s) 
    return; 
elseif ~isfield(K,'l') || isempty(K.l)
    error('Both LP and SDP cones are empty, so the problem can not be solved');
else
    K.l = K.l-1;
    if K.l == 0
        K.l = [];
    end
    K.s = 1;
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A1,c1,K1] = aggSDPcones(A,c,K0,aggConeSize)
% 
%   [A1,c1,K1] = aggSDPcones(A,c,K0,aggConeSize)
%    Input : SeDuMi format
%    aggConeSize : the block size into which small blocks are converted
%    Output: Sedumi format with larger SDP cone matrices
%

if (size(K0.s,2) > 1)
    K0.s = K0.s';
end

if size(A,1) > size(A,2)
    A = A';
end
%     if size(c,1) < size(c,2)
%         c = c';
%     end
[m,n] = size(A);
%
if ~isfield(K0,'f') || isempty(K0.f)
    fDim = 0;
else
    fDim = K0.f;
end
%
if ~isfield(K0,'l') || isempty(K0.l)
    ellDim = 0;
else
    ellDim = K0.l;
end
%
if ~isfield(K0,'q')
    qDim = 0;
else
    if size(K0.q,1) > size(K0.q,2)
        K0.q = K0.q'; % a row vector
    end
    qDim = sum(K0.q);
end
%
if ~isfield(K0,'s') || isempty(K0.s)
    K1 = K0;
    A1 = A;
    c1 = c;
    return
else
    if size(K0.s,2) > size(K0.s,1)
        K0.s = K0.s'; % a column vector
    end
    sDim = sum(K0.s .* K0.s);
end

nonSDim = fDim+ellDim+qDim;
if nonSDim == 0
    cNonSDP = [];
    AnonSDP = [];
    c0SDP = c;
    A0SDP = A;
else
    cNonSDP = c(1:nonSDim,1);
    AnonSDP = A(:,1:nonSDim);
    c0SDP = c(nonSDim+1:n,1);
    A0SDP = A(:,nonSDim+1:n);
end

K1s = [];
coneSize = 0;
for p=1:length(K0.s)
    if (K0.s(p) > aggConeSize)
        if coneSize > 0
            K1s = [K1s; coneSize; K0.s(p)];
            coneSize = 0;
        else
            K1s = [K1s; K0.s(p)];
            coneSize = 0;
        end
    elseif ((coneSize + K0.s(p) > aggConeSize))
        K1s = [K1s; coneSize];
        coneSize = K0.s(p);
    elseif ((coneSize + K0.s(p) == aggConeSize))
        K1s = [K1s; coneSize + K0.s(p)];
        coneSize = 0;
    else
        coneSize = coneSize + K0.s(p);
    end
end
if (coneSize > 0)
    K1s = [K1s; coneSize];
end

K0s = K0.s; % K0.s is a column vector
K1 = K0;
K1.s = K1s; % K1s is a column vector

tStart = tic;
% C0SDP, K0s and K1s need to be clolumn vector
[A1SDP,c1SDP] = mexAggSDPcones(A0SDP,c0SDP',K0s,K1s);
tElapsed = toc(tStart);
fprintf('Original SDP:    the max cone size, the number of cones = %3d, %3d\n',full(max(K0.s)),length(K0.s));
fprintf('Transformed SDP: the max cone size, the number of cones = %3d, %3d\n',full(max(K1.s)),length(K1.s));
% fprintf('Elapsed time for transformation = %6.2e\n',tElapsed);

c1SDP = c1SDP';

A1 = sparse([AnonSDP, A1SDP]);
c1 = sparse([cNonSDP;c1SDP]);

%    checkData(A1,b1,c1,K1);

%    sedumi(A1,b1,c1,K1);

return



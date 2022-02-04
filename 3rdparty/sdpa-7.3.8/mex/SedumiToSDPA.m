function SedumiToSDPA(filename,A,b,c,K, accuracy);
% SedumiToSDPA(filename,A,b,c,K);
%
% A converter from SeDuMi Input to SDPA sparse format
% (C) SDPA Project 2008
%
% filename : Filename for SDPA dat-s
% A,b,c,K  : SeDuMi Input 
% accuracy : printf format (e.g. %8.16e)    
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
    
    
if nargin < 6 || ~ischar(accuracy)
    accuracy = '%8.16e';
    fprintf('Accuracy is set as "%s"\n', accuracy);
end
if accuracy(1) ~= '%'
    error('accuracy must start with %% (e.g. %8.16e) \n');
end

bprint = sprintf('%s ', accuracy);
eleprint = sprintf('%%d %%d %%d %%d %s\n', accuracy);
    
% Note that primal-dual is reverse in SeDuMi
    c = -c;
    
    if isfield(K,'q') && ~isempty(K.q)
        error('Current Program cannot handle K.q');
    end

    if isfield(K,'r') && ~isempty(K.r)
        error('Current Program cannot handle K.r');
    end

    if size(b,2) ~= 1 
        b = b';
    end
    m = size(b,1);
    if size(c,2) ~= 1
        c = c';
    end
    n = size(c,1);
    if size(A,1) ~= m
        A = A';
    end
    
    fprintf('Size A[m=%d,n=%d], b[m=%d], c[n=%d] ::', ...
            size(A,1),size(A,2), m, n);
    fprintf('nnz(A) = %d, nnz(C) = %d\n', nnz(A), nnz(c));
    if size(A,1) ~= m | size(A,2) ~= n
        error('Inconsistent Size');
    end
    
    Kf = 0;
    if isfield(K,'f') && ~isempty(K.f)
        Kf = K.f;
    end
    Kl = 0;
    if isfield(K,'l') && ~isempty(K.l)
        Kl = K.l;
    end
    Ks = 0;
    if isfield(K,'s') && ~isempty(K.s)
        Ks = K.s;
        if size(Ks,2) ~= 1
            Ks = Ks';
        end
    else
        error('Cannot convert empty K.s problem');
    end
    
    fprintf('K.f = %d, K.l = %d, sum(K.s .* K.s) = %d, #K.s = %d, max(K.s) = %d\n', Kf, Kl, ...
            sum(Ks.*Ks), length(K.s), max(K.s));
    Ktotal = Kf + Kl + sum(Ks.*Ks);
    if  Ktotal ~= n
        error('Inconsistent Size K and n\n');
    end
        
    if Kf ~= 0
        fprintf(['Free Variables are divided into positive and ' ...
                 'negative part of LP cone\n']);
        Af = A(:,1:Kf);
        Al = A(:,Kf+1:Kf+Kl);
        As = A(:,Kf+Kl+1:Kf+Kl+sum(Ks.*Ks));
        Anew = [Af, -Af, Al, As];
        
        cf = c(1:Kf,1);
        cl = c(Kf+1:Kf+Kl,1);
        cs = c(Kf+Kl+1:Kf+Kl+sum(Ks.*Ks),1);
        cnew = [cf; -cf; cl; cs];
        
        Knew.f = 0;
        Knew.l = 2*Kf + Kl;
        Knew.s = K.s;
        
        A = Anew;
        c = cnew;
        K = Knew;
    end
    
    Kf = 0;
    Kl = 0;
    if isfield(K,'l') && ~isempty(K.l)
        Kl = K.l;
    end
    Ks = 0;
    if isfield(K,'s') && ~isempty(K.s)
        Ks = full(K.s);
        if size(Ks,2) ~= 1
            Ks = Ks';
        end
    end
    K.l = Kl;
    K.s = Ks;
    
    Ktotal = Kf + Kl + sum(Ks.*Ks);
    if  Ktotal ~=  size(A,2);
        error('Inconsistent Size K = %d and n = %d', Ktotal,size(A,2));
    end
    At = A';
    
    USE_MEX = 1;
    if USE_MEX
        if ~issparse(A)
            A = sparse(A);
        end
        if issparse(b)
            b = full(b);
        end
        if ~issparse(c)
            c = sparse(c);
        end
        fprintf('Writing data to %s\n', filename); 
        mexWriteSedumiToSDPA(filename, At,b,c,K, accuracy);
    else
        fid = fopen(filename,'w');
        fprintf('Writing data to %s\n', filename); 
        fprintf(fid, '%d\n',m);
        if Kl == 0
            isKl = 0;
            fprintf(fid, '%d\n',size(Ks,1));
        else
            isKl = 1;
            fprintf(fid, '%d\n',1+size(Ks,1));
            fprintf(fid, '-%d ', Kl);
        end
        fprintf(fid, '%d ', Ks);
        fprintf(fid, '\n');
        fprintf(fid, bprint , full(b));
        fprintf(fid, '\n');
        
        %  c
        
        if Kl ~= 0
            cl = c(1:Kl);
            [i,j,v] = find(cl);
            if isempty(i) ~=1
                kdummy = 0 * ones(size(i,1),1);
                ldummy = 1 * ones(size(i,1),1);
                ge = [kdummy, ldummy, i, i, v]';
                fprintf(fid, eleprint ,ge);
            end
        end
        index = Kl;
        for l=1:size(Ks,1)
            cs = c(index+1:index+Ks(l)*Ks(l));
            CS = reshape(cs,Ks(l),Ks(l));
            [i,j,v] = find(tril(CS));
            if isempty(i) ~=1
                kdummy = 0 * ones(size(i,1),1);
                ldummy = (l+isKl) * ones(size(i,1),1);
                ge = [kdummy, ldummy, i, j, v]';
                fprintf(fid, eleprint ,ge);
            end
            index = index + Ks(l) * Ks(l);
        end
        
        % A
        for k=1:m
            ak = At(:,k)';
            if Kl ~= 0
                akl = ak(1:Kl)';
                [i,j,v] = find(akl);
                if isempty(i) ~=1
                    kdummy = k * ones(size(i,1),1);
                    ldummy = 1 * ones(size(i,1),1);
                    ge = [kdummy, ldummy, i, i, v]';
                    fprintf(fid, eleprint, ge);
                end
            end
            index = Kl;
            for l=1:size(Ks,1)
                aks = ak(index+1:index+Ks(l)*Ks(l));
                AKS = reshape(aks,Ks(l),Ks(l));
                [i,j,v] = find(tril(AKS));
                if isempty(i) ~=1
                    kdummy = k * ones(size(i,1),1);
                    ldummy = (l+isKl) * ones(size(i,1),1);
                    ge = [kdummy, ldummy, i, j, v]';
                    fprintf(fid, eleprint, ge);
                end
                index = index + Ks(l) * Ks(l);
            end
        end
    end

    fprintf('Sucessfully converted\n');
    
end

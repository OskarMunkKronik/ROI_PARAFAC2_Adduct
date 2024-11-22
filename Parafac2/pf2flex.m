function [model,time] = pf2flex(X,F,options,Load)
tic
%PF2FLEX
% X         : input data mass x elution x sample
% F         : Number of components
% options   : options.MaxIter {default 1500}
%             options.ConvCrit {default 1e-9}

nonnegmeth=1;
if (nargin == 1 && ischar(X) && strcmp(X,'init'))
    detailStr = struct('ssq',[],'stopit',[],'mu',[]);
    model = struct('Loads',[],'P',[],'H',[],'detail',detailStr,'cc',[],'history',[]);
    return
end
if nargin>2
    try
        maxit = options.MaxIter;
    catch
        maxit = 1500;
    end
    try
        convcrit = options.ConvCrit;
    catch
        convcrit = 1e-6;
    end
    if (isfield(options,'show'))
        show = options.show;
    else
        show = 50;
    end
else
    convcrit = 1e-6;
    maxit = 1500;
end


sizeX=size(X);
K = sizeX(end);
if nargin<4
    % Initialize
    for i=1:length(sizeX)
        if i~=2 % The shifting mode
            Load{i} = rand(sizeX(i),F);
        else
            Load{2} = rand(size(X,2),F,K);
%            P = rand(size(X,2),F,K);
        end
    end
    H = rand(F,F);
% tic
% model = pf2flex2(X,F,options);
% Load = model.Loads;
% H = model.H;
% toc
else
%    P = rand(size(X,2),F,K);
    H = rand(F,F);
end

if ~isfield(options,'increase')
    options.increase = 1.05;
end
if ~isfield(options,'mu_lvl')
    options.mu_lvl   = 0.2;
end
if ~isfield(options,'max_mu')
    options.max_mu   = 10*norm(X(:));
end

mu = ones(1,sizeX(end));
err = 0;
for k=1:K
    % Eq 10 of original paper
    Ek = sum(sum((X(:,:,k)-Load{1}*diag(Load{3}(k,:))*Load{2}(:,:,k)').^2));
    err = err+Ek;
    mu(k) = Ek/norm(Load{2}(:,:,k),'fro')^2*(10^-1);
end

FF = NaN(maxit,3);
errold=err*2;
it=0;
P = NaN(size(Load{2}));
doIter = true(2,1);
while all(doIter) 
    it=it+1;
    errold=err;
    OldLoad=Load;
    
    % Update mu
    mu = min(mu*options.increase,options.max_mu);
    
    % Update P
    for k=1:K
        [U,~,V]=svd(H*Load{2}(:,:,k)','econ');
        P(:,:,k) = V(:,1:F)*U(:,1:F)';
    end
    
    % Update H/B*
    H = mu(1)*P(:,:,1)'*Load{2}(:,:,1);
    for k=2:K
        H = H+mu(k)*P(:,:,k)'*Load{2}(:,:,k);
    end
    H = H*diag(1./sum(H.^2));
    
    % Update A
    L = Load{2}(:,:,1)*diag(Load{end}(1,:));
    XY=X(:,:,1)*L;
    XX = L'*L;
    for k=2:K
        L = Load{2}(:,:,k)*diag(Load{end}(k,:));
        XY=XY+X(:,:,k)*L;
        XX = XX+L'*L;
    end
    if nonnegmeth==1
        pxx = pinv(XX);
        for j=1:size(XY,1)
            Aj = pxx*XY(j,:)';
            Aj(Aj<0)=0;
            A(j,:) = Aj';
        end
    elseif nonnegmeth==2
        for j=1:size(XY,1)
            A(j,:) = fnnls_old(XX,XY(j,:)')';
        end
    end
    for fa = 1:size(A,2)
        A(:,fa)=A(:,fa)/norm(A(:,fa));
    end
    Load{1}=A;
    if any(max(abs(A))<eps*1000) % Zero column
            Load{1}=A*.5+OldLoad{1}*.5;
            %disp('Correcting for zero columns in A')
    end
    if any(isnan(sum(Load{1})))
        id = find(isnan(Load{1}(:)));
        Load{1}(id)=rand(size(id));
    end

    
    % Update Bk
    Bhere = zeros(size(P,1),size(P,2));
    for k=1:K
        Ah = Load{1}*diag(Load{end}(k,:));
        Ah = [Ah ;eye(F)];
        XY = [X(:,:,k); mu(k)*H'*P(:,:,k)']'*Ah;
        XX = Ah'*Ah;
        if nonnegmeth==1
            pxx = pinv(XX);
            for j=1:size(XY,1)
                Bj = pxx*XY(j,:)';
                Bj(Bj<0)=0;
                Bhere(j,:) = Bj;
            end
            
        elseif nonnegmeth==2
            for j=1:size(XY,1)
                Bhere(j,:) = fnnls_old(XX,XY(j,:)');
            end
        end
        for fa = 1:size(Bhere,2)
            if norm(Bhere(:,fa))>eps*1000
                Bhere(:,fa)=Bhere(:,fa)/norm(Bhere(:,fa));
            end
        end
        Load{2}(:,:,k)=Bhere;
        if any(max(abs(Bhere))<eps*1000) % Zero column
            try
            Load{2}(:,:,k)=Bhere*.1+OldLoad{2}(:,:,k)*.9;
            %disp('Correcting for zero columns in B')
            catch
                3
            end
        end
    end
    
    % Update C
    for k=1:K
        AB = krb(Load{2}(:,:,k),Load{1});
        XY = vec(X(:,:,k))'*AB;
        XX = AB'*AB;
        if nonnegmeth==1
            Cj = pinv(XX)*XY';
            Cj(Cj<0)=0;
            Load{3}(k,:) = Cj;
        elseif nonnegmeth==2
            Load{3}(k,:) = fnnls_old(XX,XY');
        end    
    end
    if any(max(abs(Load{3}))<eps*1000) % Zero column
            Load{3}=Load{3}*.1+OldLoad{3}*.9;
            %disp('Correcting for zero columns in C')
    end
    
    err = 0;
    for k=1:K
        Ek=sum(sum((X(:,:,k)-Load{1}*diag(Load{3}(k,:))*Load{2}(:,:,k)').^2));
        err = err+Ek;
        if it==1 % Eq. 10 in original paper
            scalit = mu(k)*norm(Load{2}(:,:,k) - P(:,:,k)*H,'fro')^2;
            mu(k) = options.mu_lvl*Ek/scalit;
        end
    end
    if isfinite(show) && show > 0 && rem(it,show)==0 
        disp([' Fit ',num2str(err),' after ',num2str(it),' iterations'])
    end
    
    FF(it,:) = [min(mu) max(mu) err];
    
    doIter(1) = abs(err-errold)/err > convcrit;
    doIter(2) = it < maxit;
    
end
if (isfinite(show))
    if (doIter(1) && ~doIter(2)), fprintf('The algorithm did not converge after %3i iterations (fit: %8.6f; delta fit: %1.2e)\n',it,err,abs(err-errold)/err);
    else, fprintf('The algorithm converged after %3i iterations (fit: %8.6f)\n',it,err);
    end
end
detailStr = struct('ssq',err,'stopit',it,'mu',mu);
model     = struct('Loads',{Load},'P',P,'H',H,'detail',detailStr,'cc',[],'history',FF(1:it,:));
    time=toc;
if (exist('corcond','file'))
    
    dimX   = size(X);
    dimX(end) = F;
    Xtilde = NaN(dimX);
    for (k = 1:size(X,ndims(X))); Xtilde = X(:,:,k) * P(:,:,k); end
    model.cc = corcond(Xtilde,cat(2,Load(1),H,Load(end)),[],0);
    
else, model.cc = [];

end
% doplot=0;
% if doplot
%     subplot(2,2,1)
%     plot(FF(1:2,:)')
%     title('mu')
%     subplot(2,2,3)
%     plot(log(FF(3,:))')
%     title('Fits')
%     subplot(2,2,2)
%     for k=1:K
%         plot(Load{2}(:,:,k)*diag(Load{3}(k,:)))
%         hold on
%     end
%     hold off
%     shg
% end


function [x,w] = fnnls_old(XtX,Xty,tol)
%FNNLS	Non-negative least-squares.
%
% 	Adapted from NNLS of Mathworks, Inc.
%
%	x = fnnls(XtX,Xty) returns the vector X that solves x = pinv(XtX)*Xty
%	in a least squares sense, subject to x >= 0.
%	Differently stated it solves the problem min ||y - Xx|| if
%	XtX = X'*X and Xty = X'*y.
%
%	A default tolerance of TOL = MAX(SIZE(XtX)) * NORM(XtX,1) * EPS
%	is used for deciding when elements of x are less than zero.
%	This can be overridden with x = fnnls(XtX,Xty,TOL).
%
%	[x,w] = fnnls(XtX,Xty) also returns dual vector w where
%	w(i) < 0 where x(i) = 0 and w(i) = 0 where x(i) > 0.
%
%	See also NNLS and FNNLSb

%	L. Shure 5-8-87
%	Revised, 12-15-88,8-31-89 LS.
%	(Partly) Copyright (c) 1984-94 by The MathWorks, Inc.

%	Modified by R. Bro 5-7-96 according to
%       Bro R., de Jong S., Journal of Chemometrics, 1997, xx
% 	Corresponds to the FNNLSa algorithm in the paper
%
%
%	Rasmus bro
%	Chemometrics Group, Food Technology
%	Dept. Dairy and Food Science
%	Royal Vet. & Agricultural
%	DK-1958 Frederiksberg C
%	Denmark
%	rb@kvl.dk
%	http://newton.foodsci.kvl.dk/rasmus.html


%  Reference:
%  Lawson and Hanson, "Solving Least Squares Problems", Prentice-Hall, 1974.

warning off MATLAB:nearlySingularMatrix

% initialize variables
if nargin < 3
    tol = 10*eps*norm(XtX,1)*size(XtX,1);
end
[~,n] = size(XtX);
P = zeros(1,n);
Z = 1:n;
x = P';
ZZ=Z;
w = Xty-XtX*x;

% set up iteration criterion
iter = 0;
itmax = 30*n;

% outer loop to put variables into set to hold positive coefficients
while any(Z) && any(w(ZZ) > tol)
    [~,t] = max(w(ZZ));
    t = ZZ(t);
    P(1,t) = t;
    Z(t) = 0;
    PP = find(P);
    ZZ = find(Z);
    nzz = size(ZZ);
    z(PP')=(Xty(PP)'/XtX(PP,PP)');
    z(ZZ) = zeros(nzz(2),nzz(1))';
    z=z(:);
    % inner loop to remove elements from the positive set which no longer belong
    
    while any((z(PP) <= tol)) & iter < itmax
        
        iter = iter + 1;
        QQ = find((z <= tol) & P');
        alpha = min(x(QQ)./(x(QQ) - z(QQ)));
        x = x + alpha*(z - x);
        ij = find(abs(x) < tol & P' ~= 0);
        Z(ij)=ij';
        P(ij)=zeros(1,length(ij));
        PP = find(P);
        ZZ = find(Z);
        nzz = size(ZZ);
        z(PP)=(Xty(PP)'/XtX(PP,PP)');
        z(ZZ) = zeros(nzz(2),nzz(1));
        z=z(:);
    end
    x = z;
    w = Xty-XtX*x;
end



function model = pf2flex2(X,F,options,Load)

%PF2FLEX
% X         : input data mass x elution x sample
% F         : Number of components
% options   : options.MaxIter {default 1500}
%             options.ConvCrit {default 1e-9}

nonnegmeth=1;

if nargin>2
    try
        maxit=options.MaxIter;
    catch
        maxit = 1500;
    end
    try
        convcrit = options.ConvCrit;
    catch
        convcrit = 1e-9;
    end
else
    convcrit = 1e-9;
    maxit = 1500;
end


sizeX=size(X);
K = sizeX(end);
if nargin<4
    % Initialize
    for i=1:length(sizeX)
        if i~=2 % The shifting mode
            Load{i} = rand(sizeX(i),F);
        else
            Load{2} = rand(size(X,2),F,K);
            P = rand(size(X,2),F,K);
        end
    end
else
    P = rand(size(X,2),F,K);
end
H = rand(F,F);

if ~isfield(options,'increase')
    options.increase = 1.05;
end
if ~isfield(options,'mu_lvl')
    options.mu_lvl   = 0.2;
end
if ~isfield(options,'max_mu')
    options.max_mu   = 10*norm(X(:));
end

mu = ones(1,sizeX(end));
err = 0;
for k=1:K
    % Eq 10 of original paper
    Ek = sum(sum((X(:,:,k)-Load{1}*diag(Load{3}(k,:))*Load{2}(:,:,k)').^2));
    err = err+Ek;
    mu(k) = Ek/norm(Load{2}(:,:,k),'fro')^2*(10^-1);
end

FF = [];
errold=err*2;
it=0;
while (abs(err-errold)/err)>convcrit & it < maxit;
    it=it+1;
    errold=err;
    OldLoad=Load;
    
    % Update mu
    mu = min(mu*options.increase,options.max_mu);
    
    % Update P
    for k=1:K
        [U,~,V]=svd(H*Load{2}(:,:,k)','econ');
        P(:,:,k) = V(:,1:F)*U(:,1:F)';
    end
    
    % Update H/B*
    H = mu(1)*P(:,:,1)'*Load{2}(:,:,1);
    for k=2:K
        H = H+mu(k)*P(:,:,k)'*Load{2}(:,:,k);
    end
    H = H*diag(1./sum(H.^2));
    
    % Update A
    L = Load{2}(:,:,1)*diag(Load{end}(1,:));
    XY=X(:,:,1)*L;
    XX = L'*L;
    for k=2:K
        L = Load{2}(:,:,k)*diag(Load{end}(k,:));
        XY=XY+X(:,:,k)*L;
        XX = XX+L'*L;
    end
    if nonnegmeth==1
        pxx = pinv(XX);
        for j=1:size(XY,1)
            Aj = pxx*XY(j,:)';
            Aj(Aj<0)=0;
            A(j,:) = Aj';
        end
    elseif nonnegmeth==2
        for j=1:size(XY,1)
            A(j,:) = fnnls_old(XX,XY(j,:)')';
        end
    end
    for fa = 1:size(A,2)
        A(:,fa)=A(:,fa)/norm(A(:,fa));
    end
    Load{1}=A;
    if any(max(abs(A))<eps*1000) % Zero column
            Load{1}=A*.5+OldLoad{1}*.5;
            %disp('Correcting for zero columns in A')
    end
    
    % Update Bk
    Bhere = repmat(0,size(P,1),size(P,2));
    for k=1:K
        Ah = Load{1}*diag(Load{end}(k,:));
        Ah = [Ah ;eye(F)];
        XY = [X(:,:,k); mu(k)*H'*P(:,:,k)']'*Ah;
        XX = Ah'*Ah;
        if nonnegmeth==1
            pxx = pinv(XX);
            
            for j=1:size(XY,1)
                Bj = pxx*XY(j,:)';
                Bj(Bj<0)=0;
                Bhere(j,:) = Bj;
            end
            
        elseif nonnegmeth==2
            for j=1:size(XY,1)
                Bhere(j,:) = fnnls_old(XX,XY(j,:)');
            end
        end
        for fa = 1:size(Bhere,2)
            if norm(Bhere(:,fa))>eps*1000
                Bhere(:,fa)=Bhere(:,fa)/norm(Bhere(:,fa));
            end
        end
        Load{2}(:,:,k)=Bhere;
        if any(max(abs(Bhere))<eps*1000) % Zero column
            try
            Load{2}(:,:,k)=Bhere*.1+OldLoad{2}(:,:,k)*.9;
            %disp('Correcting for zero columns in B')
            catch
                3
            end
        end
    end
    
    % Update C
    for k=1:K
        AB = krb(Load{2}(:,:,k),Load{1});
        XY = vec(X(:,:,k))'*AB;
        XX = AB'*AB;
        if nonnegmeth==1
            Cj = pinv(XX)*XY';
            Cj(Cj<0)=0;
            Load{3}(k,:) = Cj;
        elseif nonnegmeth==2
            Load{3}(k,:) = fnnls_old(XX,XY');
        end    
    end
    if any(max(abs(Load{3}))<eps*1000) % Zero column
            Load{3}=Load{3}*.1+OldLoad{3}*.9;
            %disp('Correcting for zero columns in C')
    end
    
    err = 0;
    for k=1:K
        Ek=sum(sum((X(:,:,k)-Load{1}*diag(Load{3}(k,:))*Load{2}(:,:,k)').^2));
        err = err+Ek;
        if it==1 % Eq. 10 in original paper
            scalit = mu(k)*norm(Load{2}(:,:,k) - P(:,:,k)*H,'fro')^2;
            mu(k) = options.mu_lvl*Ek/scalit;
        end
    end
    if rem(it,50)==0
        disp([' Fit ',num2str(err),' after ',num2str(it),' iterations'])
    end
    
    FF=[FF [min(mu) max(mu) err]'];
end

model.Loads = Load;
model.detail.ssq = err;
model.detail.stopit = it;
model.detail.mu = mu;
model.H = H;
doplot=0;
if doplot
    subplot(2,2,1)
    plot(FF(1:2,:)')
    title('mu')
    subplot(2,2,3)
    plot(log(FF(3,:))')
    title('Fits')
    subplot(2,2,2)
    for k=1:K
        plot(Load{2}(:,:,k)*diag(Load{3}(k,:)))
        hold on
    end
    hold off
    shg
end
function AB = krb(A,B);
%KRB Khatri-Rao product
%
% The columnwise Khatri-Rao-Bro product (Harshman, J.Chemom., 2002, 198-205)
% For two matrices with similar column dimension the khatri-Rao-Bro product
% is krb(A,B) = [kron(A(:,1),B(:,1)) .... kron(A(:,F),B(:,F))]
% 
% I/O AB = krb(A,B);
%

% Copyright (C) 1995-2006  Rasmus Bro & Claus Andersson
% Copenhagen University, DK-1958 Frederiksberg, Denmark, rb@life.ku.dk
%
% This program is free software; you can redistribute it and/or modify it under 
% the terms of the GNU General Public License as published by the Free Software 
% Foundation; either version 2 of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT 
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with 
% this program; if not, write to the Free Software Foundation, Inc., 51 Franklin 
% Street, Fifth Floor, Boston, MA  02110-1301, USA.

% $ Version 1.02 $ Date 28. July 1998 $ Not compiled $
% $ Version 2.00 $ May 2001 $ Changed to array notation $ RB $ Not compiled $
% $ Version 2.01 $ May 2001 $ Error in helpfile - A and B reversed $ RB $ Not compiled $

[I,F]=size(A);
[J,F1]=size(B);

if F~=F1
   error(' Error in krb.m - The matrices must have the same number of columns')
end

AB=zeros(I*J,F);
for f=1:F
   ab=B(:,f)*A(:,f).';
   AB(:,f)=ab(:);
end


function Y = vec(X);

Y = X(:);
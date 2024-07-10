function [x,residual,iter]=GMRES_mdr(A,b,tol,maxit,x0,m,k,gs_type)
%
%  x = approximate solution to A*x=b
%
%  residual = history of l2 residual norm for each restart cycle
%
%  iter = iteration counts
%
%  tol = desired accuracy of the solution measured
%  by relative residual error
%  If tol is [] then we take tol = 1e-6.
%
%  maxit = the maximum number of Arnoldi iterations
%  If maxit is [] then maxit = 200.
%
%  x0 = the initial guess.
%  If x0 is [] then x0 is defined by n x 1 zero vector.
%
%  m = restart iteration
%  If m is [] then m = 50.
%
%  k = the number of deflated vectors
%  If k is [] then k = 0.
%
%  gs_type = type of Gram-Schmidt process
%  gs_type = 0 --> CGS (by default)
%  gs_type = 1 --> CGS2
%  gs_type = 2 --> RGS+CGS
%  gs_type = 3 --> RGS+MGS
%
%  Copyright (c) 2024, Yongseok Jang. 
%
%  This program is free software: you can redistribute it and/or modify 
%  it under the terms of the GNU Lesser General Public License as 
%  published by the Free Software Foundation, either version 3 of 
%  the License, or (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful, 
%  but WITHOUT ANY WARRANTY; without even the implied warranty of 
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
%  See the GNU Lesser General Public License for more details.
%
%  You should have received a copy of the GNU Lesser General Public License 
%  along with this program. If not, see <https://www.gnu.org/licenses/>. 

if (nargin < 2)
    error(message('Not enough input variables'));
end
n=length(b);
if (nargin < 3) || isempty(tol)
        tol = 1e-6;
end
if (nargin < 4) || isempty(maxit)
        maxit = 200;
end
if (nargin < 5) || isempty(x0)
        x0 = zeros(n,1);
end
if (nargin < 6) || isempty(m)
        m = 50;
end
if (nargin < 7) || isempty(k)
        k = 0;
end
if (nargin < 8) || isempty(gs_type)
        gs_type = 0;
end

% convert a matrix A to a function
if isa(A,'sparse')||isa(A,'double')
    A = @(x) A*x;
elseif ~isa(A,'function_handle')
    error(message('randgmres: l.h.s. must be a square matrix or function handle'));
end

if gs_type > 1
    t = min(n, ceil(2*m*log(n)/log(m))); % Sketch size
    Theta = SRHT(n,t);
else
    Theta = [];
end
r0=b-A(x0);

beta=norm(r0);
nb=beta;
c=[beta,zeros(1,m)]';

iter=zeros(maxit,1);
residual=zeros(maxit,1);
v0=r0/beta;
residual(1)=norm(r0);
iter(1)=0;

% first m-Arnoldi iterations
[V,H] = Arnoldi(A,n,m,v0,gs_type,Theta);

y=H\c;
x=x0+V(:,1:m)*y;
r0=b-A(x);

beta=norm(r0);
residual(2)=beta/nb;
iter(2)=m;
it=1;
if beta<tol*nb
    residual=residual(1:it+1);
    iter=iter(1:it+1);
    return
end

ek=zeros(m+1,1);
ek(k+1)=1;
W=V(:,1:m);
while beta>tol*nb&&iter(it+1)<maxit
    it=it+1;
    % new computation with harmonic Ritz forms
    if k~=0
        T1=H'*H;
        if it==2
            T2=eye(m,m);
        else
%             T2=[Wk'*Wk,zeros(k,m-k);zeros(m-k,k),eye(m-k,m-k)];
            T2=W'*W;
        end
        [E,D]=eig(T1,T2);
        e=diag(D);
        [~,idx]=sort(e);
        id=idx(1:k);
        P=E(:,id);
        [Q,R]=qr(H*P,0);
        invR=pinv(R);
        Vk=V*Q;
        Wk=W*P*invR;
    end
    
    % (m-k) Arnoldi iterations
    if k==0
        v0=r0/beta;
        [V,H] = Arnoldi(A,n,m,v0,gs_type,Theta);
        W=V(:,1:m);
    else
        v0=r0/beta;
        [V2,H2] = Arnoldi_mdr(A,m,gs_type,Theta,k,Vk,v0);
    end    
    if k==0
        c=[beta,zeros(1,m)]';
    else
        H=[eye(k,k),Vk'*A(V2(:,1:end-1));zeros(m-k+1,k),H2];
        V=[Vk,V2];
        W=[Wk,V2(:,1:end-1)];
        c=beta*ek;
    end
    y=H\c;
    x=x+W*y;
    r0=b-A(x);    
    beta=norm(r0);   
    residual(it+1)=norm(r0)/nb;
    if it==1
        iter(it+1)=iter(it)+(m-k)+m;
    else
        iter(it+1)=iter(it)+(m-k);
    end
end
residual=residual(1:it+1);
iter=iter(1:it+1);
end
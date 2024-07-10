function [x,residual,iter,condnum,loss]=GMRES_dr(A,b,tol,maxit,x0,m,k,gs_type,ilu0)
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
%  ilu0 = switch applying ILU0 preconditioning
%  0: off,  1: on
%  If ilu0 is [] then ilu0 = 0.
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
if (nargin < 9) || isempty(ilu0)
        ilu0 = 0;
end

% define ILU preconditioner
if ilu0 == 1
    [L,U] = ilu(A);
    prec = @(x) U\(L\x);
else
    prec = @(x) x;
end

% convert a matrix A to a function
if isa(A,'sparse')||isa(A,'double')
    A0 = @(x) A*x; %A0 is the original system
elseif ~isa(A,'function_handle')
    error(message('randgmres: l.h.s. must be a square matrix or function handle'));
end

A = @(x) A0(prec(x));

if gs_type > 1
    t = min(n, ceil(2*m*log(n)/log(m))); % Sketch size
    fprintf('dimension of sketch = %d\n',t);
    Theta = SRHT(n,t);
else
    Theta = [];
end
r0=b-A0(x0);

beta=norm(r0);
nb=beta;
c=[beta,zeros(1,m)]';

iter=zeros(maxit,1);
residual=zeros(maxit,1);
condnum=zeros(maxit,1);
loss=zeros(maxit,1);

v0=r0/beta;
residual(1)=norm(r0);
iter(1)=0;
condnum(1)=1;
loss(1)=0;

% first m-Arnoldi iterations
[V,H] = Arnoldi(A,n,m,v0,gs_type,Theta);
y=H\c;
x=x0+prec(V(:,1:m)*y);
r0=b-A0(x);

beta=norm(r0);
residual(2)=beta/nb;
iter(2)=m;

condnum(2)=cond(V);
loss(2)=norm(eye(m+1,m+1)-V'*V);

it=1;
if beta<tol*nb
    residual=residual(1:it+1);
    iter=iter(1:it+1);
    return
end

em=zeros(m,1);
em(m)=1;


while beta>tol*nb&&iter(it+1)<maxit
    it=it+1;
    % new computation with harmonic Ritz forms
    if k~=0
        rho=c-H*y;
        H2=H(1:m,1:m);
        h=H(m+1,m);
        f=pinv(H2')*em;
        T=H2+(h*h*f*(em'));
        [E,D]=eig(T);
        e=diag(D);
        [~,idx]=sort(abs(e),"ascend");
        id=idx(1:k);
%         disp('Harmonic Ritz values')
%         disp(abs(e(id)))
%         disp('---------------')
%         pause
        G=E(:,id);
        P2=[G;zeros(1,k)];
        AP=[P2,rho];
        
        [QQ,~]=qr(AP,0);
    
        V_n=V*QQ;
        H_n=QQ'*H*QQ(1:m,1:k);
        V(:,1:k+1)=V_n;
        H(1:k+1,1:k)=H_n;
        
    end
    
    % (m-k) Arnoldi iterations
    if k==0
        v0=r0/beta;
        [V,H] = Arnoldi(A,n,m,v0,gs_type,Theta);
    else
        [V,H] = Arnoldi_dr(A,m,gs_type,Theta,k,V,H);
    end    
    if k==0
       c=[beta,zeros(1,m)]';
    else
       c=[QQ'*rho;zeros(m-k,1)];
    end
    y=H\c;
    x=x+prec(V(:,1:m)*y);
    r0=b-A0(x);    
    beta=norm(r0);   
    residual(it+1)=norm(r0)/nb;


    condnum(it+1)=cond(V);
    loss(it+1)=norm(eye(m+1,m+1)-V'*V);

    if it==1
        iter(it+1)=iter(it)+(m-k)+m;
    else
        iter(it+1)=iter(it)+(m-k);
    end
end
residual=residual(1:it+1);
iter=iter(1:it+1);

condnum=condnum(1:it+1);
loss=loss(1:it+1);

end
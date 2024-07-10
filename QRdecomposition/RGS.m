function [Q,R]=RGS(W,k)
[n,m]=size(W);
Q=zeros(n,m);

R=zeros(m,m);
S=zeros(k,m);

Theta = SRHT(n,k);

% orthogonalization with RGS
for i=1:m
    q=W(:,i);
    p=Theta(q);
    R(1:i-1,i)=S(:,1:i-1)\p;

%     [R(1:i-1,i),~]=lsqr(S(:,1:i-1),p,1e-16);
    for j=1:i-1
        q=q-Q(:,j)*R(j,i);
    end
    s=Theta(q);
    R(i,i)=norm(s);
    S(:,i)=s/R(i,i);
    Q(:,i)=q/R(i,i);
end

end
function [Q,R]=RGS2_CGS2(W,k)
[n,m]=size(W);
Q=zeros(n,m);

R1=zeros(m,m);
R2=zeros(m,m);
R=zeros(m,m);
S=zeros(k,m);

Theta = SRHT(n,k);

for i=1:m
    u=W(:,i);
    p=Theta(u);
    R1(1:i-1,i)=(S(:,1:i-1)\p);
    
    for j=1:i-1
        u=u-Q(:,j)*R1(j,i);
    end

    v=u;
    for j=1:i-1
        R2(j,i)=Q(:,j)'*v;
    end
    for j=1:i-1
        v=v-Q(:,j)*R2(j,i);
    end
    R(i,i)=norm(v);
    Q(:,i)=v/R(i,i);
    S(:,i)=Theta(Q(:,i));
    R(1:i-1,i)=R1(1:i-1,i)+R2(1:i-1,i);
end

end
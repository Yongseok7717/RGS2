function [Q,R]=MGS2(W)
[n,m]=size(W);

Q=zeros(n,m);
R=zeros(m,m);
R1=zeros(m,m);
R2=zeros(m,m);

% orthogonalization with MGS2
for j=1:m
    v=W(:,j);

    for i=1:j-1
        R1(i,j)=Q(:,i)'*v;
        v=v-R1(i,j)*Q(:,i);
    end

    w=v;
    for i=1:j-1
        R2(i,j)=Q(:,i)'*w;
        w=w-R2(i,j)*Q(:,i);
    end
    R(1:j-1,j)=R1(1:j-1,j)+R2(1:j-1,j);
    R(j,j)=norm(w);
    Q(:,j)=w/R(j,j);
end

end
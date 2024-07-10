function [Q,R]=CGS(W)
[n,m]=size(W);

Q=zeros(n,m);
R=zeros(m,m);

% orthogonalization with CGS
for j=1:m
    v=W(:,j);

    for i=1:j-1
        R(i,j)=Q(:,i)'*v;
    end
    for i=1:j-1
        v=v-R(i,j)*Q(:,i);
    end
    R(j,j)=norm(v);
    Q(:,j)=v/R(j,j);
end

end
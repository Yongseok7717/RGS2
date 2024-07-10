function [Q,R]=MGS(W)
[n,m]=size(W);

Q=zeros(n,m);
R=zeros(m,m);

% orthogonalization with MGS
for j=1:m
    v=W(:,j); 
    for i=1:j-1
        R(i,j)=Q(:,i)'*v;
        v=v-R(i,j)*Q(:,i);
    end
    R(j,j)=norm(v);
    Q(:,j)=v/R(j,j);
end

end
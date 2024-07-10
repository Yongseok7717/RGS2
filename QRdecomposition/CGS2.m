function [Q,R]=CGS2(W)
[n,m]=size(W);

Q=zeros(n,m);
R=zeros(m,m);
R1=zeros(m,m);
R2=zeros(m,m);

% orthogonalization with CGS2
for j=1:m
    v=W(:,j);

    for i=1:j-1
        R1(i,j)=Q(:,i)'*v;
    end
    for i=1:j-1
        v=v-R1(i,j)*Q(:,i);
    end

    w=v;
    for i=1:j-1
        R2(i,j)=Q(:,i)'*w;
    end
    for i=1:j-1
        w=w-R2(i,j)*Q(:,i);
    end
    
    R(1:j-1,j)=R1(1:j-1,j)+R2(1:j-1,j);
    R(j,j)=norm(w);
    Q(:,j)=w/R(j,j);
    if cond(Q(:,j))>1
        disp(j)
        break
    end
end

end

function [V,H] = Arnoldi(A,n,m,v0,gs_type,Theta)

V=zeros(n,m+1);
H=zeros(m+1,m);

V(:,1)=v0;

if gs_type == 0 % CGS
    for j=1:m
        w=A(V(:,j));
        for i=1:j
            H(i,j)=V(:,i)'*w;
        end
        for i=1:j
            w=w-H(i,j)*V(:,i);
        end
        
        H(j+1,j)=norm(w);
        
        if H(j+1,j)==0
            V(:,j+1)=0;
        else
            V(:,j+1)=w/H(j+1,j);
        end
    end
elseif gs_type == 1 % CGS2
    tempH=H;
    tempH2=H;
    for j=1:m
        w=A(V(:,j));
        for i=1:j
            tempH(i,j)=V(:,i)'*w;
        end
        for i=1:j
            w=w-tempH(i,j)*V(:,i);
        end
        ww=w;
        for i=1:j
            tempH2(i,j)=V(:,i)'*ww;
        end
        for i=1:j
            ww=ww-tempH2(i,j)*V(:,i);
        end
        H(1:j,j)=tempH(1:j,j)+tempH2(1:j,j);
        H(j+1,j)=norm(ww);
        
        if H(j+1,j)==0
            V(:,j+1)=0;
        else
            V(:,j+1)=ww/H(j+1,j);
        end
    end
elseif gs_type == 2 % RGS+CGS
    w=v0;
    p=Theta(w);
    t=length(p);
    S=zeros(t,m+1);
    S(:,1)=p;

    tempH=H;
    tempH2=H;
    for j=1:m
        w=A(V(:,j));
        p=Theta(w);
        tempH(1:j,j)=S(:,1:j)\p;
        
        for i=1:j
            w=w-V(:,i)*tempH(i,j);
        end

        ww=w;
        for i=1:j
            tempH2(i,j)=V(:,i)'*ww;
        end
        for i=1:j
            ww=ww-tempH2(i,j)*V(:,i);
        end
        H(1:j,j)=tempH(1:j,j)+tempH2(1:j,j);
        H(j+1,j)=norm(ww);
        
        if H(j+1,j)==0
            V(:,j+1)=0;
        else
            V(:,j+1)=ww/H(j+1,j);
        end
        S(:,j+1)=Theta(V(:,j+1));
    end
else % RGS+MGS
    w=v0;
    p=Theta(w);
    t=length(p);
    S=zeros(t,m+1);
    S(:,1)=p;

    tempH=H;
    tempH2=H;
    for j=1:m
        w=A(V(:,j));
        p=Theta(w);
        tempH(1:j,j)=S(:,1:j)\p;
        
        for i=1:j
            w=w-V(:,i)*tempH(i,j);
        end

        ww=w;
        for i=1:j
            tempH2(i,j)=V(:,i)'*ww;
            ww=ww-tempH2(i,j)*V(:,i);
        end
        H(1:j,j)=tempH(1:j,j)+tempH2(1:j,j);
        H(j+1,j)=norm(ww);
        
        if H(j+1,j)==0
            V(:,j+1)=0;
        else
            V(:,j+1)=ww/H(j+1,j);
        end
        S(:,j+1)=Theta(V(:,j+1));
    end
end
function [condnum,loss] = computeStab(Q)
 [~,m]=size(Q);
 sz=m/50;
 condnum=zeros(sz,1);
 loss=zeros(sz,1);
 for i=1:sz
     condnum(i)=cond(Q(:,1:i*50));
     loss(i)=norm(eye(i*50,i*50)-Q(:,1:i*50)'*Q(:,1:i*50));
 end
end

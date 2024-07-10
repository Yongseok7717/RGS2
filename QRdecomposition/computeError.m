function [error] = computeError(Q,R,W)
 [~,m]=size(Q);
 sz=m/50;
 error=zeros(sz,1);
 for i=1:sz
     error(i)=norm(W(:,1:i*50)-Q(:,1:i*50)*R(1:i*50,1:i*50));
 end
end
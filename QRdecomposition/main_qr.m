clear
clc
f=@(x,y) sin(10*(x+y))./(cos(100*(y-x))+1.1);
n=1e6;
m=500;
t = min(n, ceil(2*m*log(n)/log(m)));
fprintf('Sketch size = %d\n',t)
x=linspace(0,1,n)';
y=linspace(0,1,m)';
[X,Y]=meshgrid(x,y);

W=f(X,Y)';


% Perform QR decomposition on W w.r.t GS process
tic
[Q_rgs,R_rgs]=RGS(W,t);
timeRGS=toc;
fprintf('Elapsed time for RGS = %.5f\n',timeRGS)

tic
[Q_cgs,R_cgs]=CGS(W);
timeCGS=toc;
fprintf('Elapsed time for CGS = %.5f\n',timeCGS)

tic
[Q_mgs,R_mgs]=MGS(W);
timeMGS=toc;
fprintf('Elapsed time for MGS = %.5f\n',timeMGS)

tic
[Q_rgs_cgs,R_rgs_cgs]=RGS2_CGS2(W,t);
timeRGSCGS=toc;
fprintf('Elapsed time for RGS+CGS = %.5f\n',timeRGSCGS)

tic
[Q_rgs_mgs,R_rgs_mgs]=RGS2_MGS2(W,t);
timeRGSMGS=toc;
fprintf('Elapsed time for RGS+MGS = %.5f\n',timeRGSMGS)

tic
[Q_cgs2,R_cgs2]=CGS2(W);
timeCGS2=toc;
fprintf('Elapsed time for CGS2 = %.5f\n',timeCGS2)

tic
[Q_mgs2,R_mgs2]=MGS2(W);
timeMGS2=toc;
fprintf('Elapsed time for MGS2 = %.5f\n',timeMGS2)

%% Compute condition numbers
[cond_rgs,~]=computeStab(Q_rgs);
[cond_cgs,loss_cgs]=computeStab(Q_cgs);
[cond_cgs2,loss_cgs2]=computeStab(Q_cgs2);
[cond_mgs,loss_mgs]=computeStab(Q_mgs);
[cond_mgs2,loss_mgs2]=computeStab(Q_mgs2);
[cond_rgscgs,loss_rgscgs]=computeStab(Q_rgs_cgs);
[cond_rgsmgs,loss_rgsmgs]=computeStab(Q_rgs_mgs);

% Compute approximation errors
normW=zeros(m/50,1);
condW=zeros(m/50,1);
for i=1:m/50
    normW(i)=norm(W(:,1:50*i));
    condW(i)=cond(W(:,1:50*i));
end

[error_rgs]=computeError(Q_rgs,R_rgs,W)./normW;
[error_cgs]=computeError(Q_cgs,R_cgs,W)./normW;
[error_cgs2]=computeError(Q_cgs2,R_cgs2,W)./normW;
[error_mgs]=computeError(Q_mgs,R_mgs,W)./normW;
[error_mgs2]=computeError(Q_mgs2,R_mgs2,W)./normW;
[error_rgscgs]=computeError(Q_rgs_cgs,R_rgs_cgs,W)./normW;
[error_rgsmgs]=computeError(Q_rgs_mgs,R_rgs_mgs,W)./normW;


%% Graphics
% Stability figure
iter=1:m/50;
iter=50*iter;

fig1=figure;
semilogy(iter,cond_cgs,'-*',iter,cond_cgs2,'-s',...
    iter,cond_mgs,'-o',iter,cond_mgs2,'-d',...
    iter,cond_rgs,iter,cond_rgscgs,'-<',iter,cond_rgsmgs,'-+');
legend('CGS','CGS2','MGS','MGS2',...
    'RGS','RGS2C','RGS2M','Location','best')
xlim([50,m])
ylim([0,1e5])
title('Condition numbers of $$Q_i$$','interpreter','latex' )
xlabel('$$i$$','interpreter','latex')

saveas(fig1,'cond_num.png')

fig2=figure;
semilogy(iter,loss_cgs,'-*',iter,loss_cgs2,'-s',...
    iter,loss_mgs,'-o',iter,loss_mgs2,'-d',...
    iter,loss_rgscgs,'-<',iter,loss_rgsmgs,'-+');
legend('CGS','CGS2','MGS','MGS2',...
    'RGS2C','RGS2M','Location','best')
xlim([50,m])
ylim([0,1e3])
title('Loss of orthogonality $$||I_i-(Q_i)^TQ_i||$$','interpreter','latex' )
xlabel('$$i$$','interpreter','latex')

saveas(fig2,'loss_ortho.png')

% Approximation error

fig3=figure;
semilogy(iter,error_cgs,'-*',iter,error_cgs2,'-s',...
    iter,error_mgs,'-o',iter,error_mgs2,'-d',...
    iter,error_rgs,iter,error_rgscgs,'-<',iter,error_rgsmgs,'-+');
legend('CGS','CGS2','MGS','MGS2',...
    'RGS','RGS2C','RGS2M','Location','best')
xlim([50,m])
title('Error $$|| W_i-Q_iR_i||/|| W_i||$$','interpreter','latex')
xlabel('$$i$$','interpreter','latex')

saveas(fig3,'approx_err.png')


fig4=figure;
semilogy(iter,condW);
title('Condition numbers of $$W_i$$','interpreter','latex' )
xlabel('$$i$$','interpreter','latex')

saveas(fig1,'cond_numW.png')

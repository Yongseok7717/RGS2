clear
clc

load('ML_Geer.mat');
A=Problem.A;
b=ones(length(A),1);
b=A*b;
b=b/norm(b);

tol=1e-8;
nRestart=400;
nCycle=10;
maxit=nRestart*nCycle;

nDeflation=40;

%%CGS based
tic
[~,res,it,condNum,loss]=GMRES_dr(A,b,tol,maxit,[],nRestart,0,0,1);
timeCGS=toc;
fprintf('Elapsed time for CGS GMRES = %.5f\n',timeCGS)

tic
[~,res_dr,it_dr,condNumdr,lossdr]=GMRES_dr(A,b,tol,maxit,[],nRestart,nDeflation,0,1);
timeCGSdr=toc;
fprintf('Elapsed time for CGS GMRESDR = %.5f\n',timeCGSdr)


%%RGS2C based
tic
[~,resc,itc,condNumc,lossc]=GMRES_dr(A,b,tol,maxit,[],nRestart,0,2,1);
timeRGS2C=toc;
fprintf('Elapsed time for RGS2C GMRES = %.5f\n',timeRGS2C)

tic
[~,res_drc,it_drc,condNumdrc,lossdrc]=GMRES_dr(A,b,tol,maxit,[],nRestart,nDeflation,2,1);
timeRGS2Cdr=toc;
fprintf('Elapsed time for RGS2C GMRESDR = %.5f\n',timeRGS2Cdr)


%%RGS2C based
tic
[~,resm,itm,condNumm,lossm]=GMRES_dr(A,b,tol,maxit,[],nRestart,0,3,1);
timeRGS2M=toc;
fprintf('Elapsed time for RGS2M GMRES = %.5f\n',timeRGS2M)

tic
[~,res_drm,it_drm,condNumdrm,lossdrm]=GMRES_dr(A,b,tol,maxit,[],nRestart,nDeflation,3,1);
timeRGS2Mdr=toc;
fprintf('Elapsed time for RGS2M GMRESDR = %.5f\n',timeRGS2Mdr)

% 
% 
% %%Graphic
fig1=figure;
semilogy(it,res,'-',it_dr,res_dr,':',...
    itc,resc,'-*',...
    it_drc,res_drc,'-s',...
    itm,resm,'-d',it_drm,res_drm,'-o')
legend('GMRES(400)','GMRES DR(400,40)',...
'RGS2C GMRES(400)','RGS2C GMRES DR(400,40)',...
'RGS2M GMRES(400)','RGS2M GMRES DR(400,40)')
xlabel('Iterations')
ylabel('Relative residual norms')

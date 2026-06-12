clearvars;
clc;
format shorte;
% Spatial step
a=0; b=1; n =128;  h=(b-a)/n;
% Trial points
[x,y]=meshgrid(a:h:b);
X=[x(:) y(:)];
% Boundary points
indb1=find(X(:,1)==a);                        Xb1=X(indb1,:); Nb1=length(Xb1);
indb2=find(X(:,1)==b);                        Xb2=X(indb2,:); Nb2=length(Xb2);
indb3=find(X(:,2)==a & X(:,1)>a & X(:,1)<b);  Xb3=X(indb3,:); Nb3=length(Xb3);
indb4=find(X(:,2)==b & X(:,1)>a & X(:,1)<b);  Xb4=X(indb4,:); Nb4=length(Xb4);
Xb=[Xb1;Xb2;Xb3;Xb4]; Nb=length(Xb);
% Interior points
indi=find(X(:,1)>a & X(:,1)<b & X(:,2)>a & X(:,2)<b);  Xi=X(indi,:);  Ni=length(Xi);
% Whole points
X=[Xb;Xi];
% Constant Parameter
nu=0.0005;
% Final time & time step
T=pi/2; MT=4000; dt=T/MT;
% Exact solution 
xbar=@(x,y,t) x.*cos(4*t)+y.*sin(4*t); ybar=@(x,y,t) -x.*sin(4*t)+y.*cos(4*t);
sigma2=0.002; x0=0.5; y0=0.75;
Uex=@(x,y,t) (sigma2./(sigma2+4*nu*t)).*exp(-((xbar(x,y,t)-x0).^2+(ybar(x,y,t)-y0).^2)./(sigma2+4*nu*t));
% Initial condition
U0=Uex(X(:,1),X(:,2),0)+eps;
% Differentiation matrices
AA=Weight_Rational_RBF_PU_2D_CD([X(:,1),X(:,2)],[X(:,1),X(:,2)],a,b,n,'0',U0);
AAxx=Weight_Rational_RBF_PU_2D_CD([X(:,1),X(:,2)],[X(:,1),X(:,2)],a,b,n,'xx+yy',U0);
Ax=Weight_Rational_RBF_PU_2D_CD([X(:,1),X(:,2)],[X(:,1),X(:,2)],a,b,n,'1x',U0);
Ay=Weight_Rational_RBF_PU_2D_CD([X(:,1),X(:,2)],[X(:,1),X(:,2)],a,b,n,'1y',U0);
% Right hand-side of equation
R=@(u,ux,uy,nabla2u,x,y) nu*(nabla2u)+(4*y).*ux-(4*x).*uy;
for nt=1:MT
    % Fourht-order Runge-Kutta method
    t=nt*dt
    K1=R(AA*U0,Ax*U0,Ay*U0,AAxx*U0,X(:,1),X(:,2));
    K2=R(AA*(U0+(dt*K1/2)),Ax*(U0+(dt*K1/2)),Ay*(U0+(dt*K1/2)),AAxx*(U0+(dt*K1/2)),X(:,1),X(:,2));
    K3=R(AA*(U0+(dt*K2/2)),Ax*(U0+(dt*K2/2)),Ay*(U0+(dt*K2/2)),AAxx*(U0+(dt*K2/2)),X(:,1),X(:,2));
    K4=R(AA*(U0+(dt*K3)),Ax*(U0+(dt*K3)),Ay*(U0+(dt*K3)),AAxx*(U0+(dt*K3)),X(:,1),X(:,2));
    U0(Nb+1:N)=U0(Nb+1:N)+(dt/6)*(K1(Nb+1:N)+2*(K2(Nb+1:N)+K3(Nb+1:N))+K4(Nb+1:N));
    U0(1:Nb)=Uex(Xb(:,1),Xb(:,2),t);
    if isnan(U0)
        break;
    end
end
% Compute relative error
ER=norm(U0-Uex(X(:,1),X(:,2),t),2)/norm(Uex(X(:,1),X(:,2),t),2)
% Draw the approximate solution at test points
sz=[100 100];   % surface grid parameters
[ll,tt]=meshgrid(linspace(a,b,sz(2)),linspace(a,b,sz(1)));
xx=[ll(:) tt(:)];
% **************** Interpolate to the grid ************* %
AI=Weight_Rational_RBF_PU_2D_MHD([X(:,1),X(:,2)],xx,a,b,n,'0',U0);
yy=reshape(xx(:,2),sz);
xx=reshape(xx(:,1),sz);
% Approximate solution 
uapp=AI*U0;
indapp=find(uapp<0);
uapp(indapp)=0;
uu=reshape(uapp,sz);
mesh(xx,yy,uu)
shading interp;
str = ['t=',num2str(t)];
title(str);
colormap jet;

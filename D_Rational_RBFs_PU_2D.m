% This program gives the approximate solutions to the 2D functions and their first derivatives via the direct 
% rational radial basis functions partition of unity (DRRBF-PU) method. Two different examples are written in this code.
% The users can add new examples (see Lines 8-23) to check the results.
% The users can utilize the quasi-uniform nodes (for example, Halton nodes) by changing Lines 29 and 30.
% To run this code, it is also necessary to consider the other function "Weight_Rational_RBF_PU_2D.m".
clearvars;
clc;
format shorte;
Tfunc='R2'; % Tfunc: is the test function, i.e., 'R1' (first example) or 'R2' (second example) given in the paper;
% Test function
switch (Tfunc)
    case ('R1')
    u=@(xx,yy) atan(125*(sqrt((xx-1.5).^2+(yy-0.25).^2)-0.92));
    syms xx yy;
    ux=inline(diff(u(xx,yy),xx),'xx','yy');
    uy=inline(diff(u(xx,yy),yy),'xx','yy');
    case ('R2')
     u=@(xx,yy) tan(9*(yy-xx)+1)./(tan(9)+1);
     syms xx yy;
     ux=inline(diff(u(xx,yy),xx),'xx','yy');
     uy=inline(diff(u(xx,yy),yy),'xx','yy');
    otherwise
        error ('No example is found, but the user can add a new example...')
end
% Mesh generation (trial points) (here is uniform nodes)
n=input([' Enter the number of points in each direction (recommended n>=32 (32, 64, 128 and 256 were used) ....' ...
    'to get more accurate results):']);
disp('Please wait to run the program...')
a=0; b=1; h=(b-a)/n;
[x,y]=meshgrid(a:h:b);
X=[x(:) y(:)]; N=length(X);
% Evaluate the approximate solution at test points
sz=[100 100];
[ll,tt]=meshgrid(linspace(a,b,sz(1)),linspace(a,b,sz(2)));
Xc=[ll(:),tt(:)];
% Differentiation matrix (call subroutine "Weight_Rational_RBF_PU_2D.m")
Order=
Dx=Weight_Rational_RBF_PU_2D([X(:,1),X(:,2)],[Xc(:,1),Xc(:,2)],a,b,n,Order,Tfunc);
% plot of approximate and exact solution
f=u(X(:,1),X(:,2));
fre=Dx*f;
uappR=reshape(fre,sz);
yy=reshape(Xc(:,2),sz);
xx=reshape(Xc(:,1),sz);
uex=reshape(ux(Xc(:,1),Xc(:,2)),sz);
% Figure 1
subplot(1,2,1),mesh(xx,yy,uappR) % Approximate
xlabel('x')
ylabel('y')
zlabel('$\frac{\partial f(x,y)}{\partial x}$')
title('Approximate solution')
% Figure 2
subplot(1,2,2),mesh(xx,yy,uex) % Exact
xlabel('x')
ylabel('y')
zlabel('$\frac{\partial f(x,y)}{\partial x}$')
title('Exact solution')
% Errors (Relative error in two-norm)
ER=norm(uex-uappR,2)/norm(uex,2)

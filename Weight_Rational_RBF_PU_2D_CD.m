function A=Weight_Rational_RBF_PU_2D_CD(X,Xe,a,b,n,Order,u)
% Inputs: X is a set of trial points; 
%         Xe is a set of test points;
%         a and b are the beginning and end of the interval (Domain is the square [a,b]^2)
%         n=(b-a)/h for choosing the number of patch (Npu=(n/4)^2)
%         Order: '0' for the approximation; 
%                '1x' for the first derivative in directin x;
%                '1y' for the first derivative in directin y;
%                'xx+yy'for the second derivative (Laplace operator);

% Output: A is the differentiation matrix;

% RBF (Matern kernel)
rbf=@(r,c)       exp(-c.*r).*(15+15*(c.*r)+6*(c.*r).^2+(c.*r).^3);
dxrbf=@(r,c,x)  -c.^2.*x.*exp(-c.*r).*((c.*r).^2+3*(c.*r)+3);
dyrbf=@(r,c,y) -c.^2.*y.*exp(-c.*r).*((c.*r).^2+3*(c.*r)+3);
d2rbf=@(r,c)    c.^2.*exp(-c.*r).*((c.*r).^3-3*(c.*r)-3)-c.^2.*exp(-c.*r).*((c.*r).^2+3*(c.*r)+3);
% weight function (the PU weights in the approximation)
wf=@(r,e)       (1-e*r).^6.*(35*(e*r).^2+18*(e*r)+3);
% Number of trial points
N=length(X); 
% Number of test points
M=length(Xe);
% Number of patches
npu=n/4; xx=linspace(a,b,npu); [xl,yl]=meshgrid(xx);
Xl=[xl(:) yl(:)];             
Nc=length(Xl);
% Specific parameters related to the proposed method
delta=(b-a)/npu;
% Shape parameter
c=35;
% radius of patch 
rl=delta; ep=1;% (used in PU weights)
% Matrix related to the PU weights
SEM=zeros(M,Nc);
for i=1:M
    ind=find(sqrt((Xe(i,1)-Xl(:,1)).^2+(Xe(i,2)-Xl(:,2)).^2)<rl);
    Y=Xl(ind,:);
    rpu=sqrt((Xe(i,1)-Y(:,1)').^2+(Xe(i,2)-Y(:,2)').^2)/rl;
    SEM(i,ind)= wf(rpu,ep)./sum(wf(rpu,ep));
end
% Matrix
A=spalloc(M,N,25*N);
% Loop over each patch
for l=1:Nc
    xxl=Xl(l,1); yyl=Xl(l,2);
    % Find trial points in patch l
    indl=find(sqrt((xxl-X(:,1)).^2+(yyl-X(:,2)).^2)<delta);
    Yl=X(indl,:); 
    Nl=length(Yl);
    % Compute the interpolation matrix 
    rrl=sqrt((Yl(:,1)-Yl(:,1)').^2+(Yl(:,2)-Yl(:,2)').^2);
    RBFI=rbf(rrl,c); % A_{phi,X_{l}}
    
    % Compute the interpolation matrix
    mu=0.00000001;
    Aphi=RBFI+mu*eye(Nl);
    f=u(indl);
    D=diag(f);

    % Cholesky factorization
    gam=0;
    L=chol(Aphi,'lower');
    Ainv=L'\(L\eye(Nl))+gam*eye(Nl);

    K=1/norm(f,2)^2;
    S=diag(1./(K*f.^2+1))*(K*(D'*(Ainv)*D+Ainv));

    opts.tol=1e-6;
    [q,~]= eigs(S,Nl,'smallestreal',opts);
    q=q(:,1);

    % Find coefficients alpha and beta
    beta=L'\(L\q);
    AphiR=(1./(Aphi*beta)).*(Aphi).*(1./(Aphi*beta))';
    
    % Find test points in patch l
    indle=find(sqrt((xxl-Xe(:,1)).^2+(yyl-Xe(:,2)).^2)<delta);
    Yle=Xe(indle,:); 
    % Compute the interpolation matrix and its rational approximation at test points
    rle=sqrt((Yle(:,1)-Yl(:,1)').^2+(Yle(:,2)-Yl(:,2)').^2);
    RBFe=rbf(rle,c);
    Aphic=RBFe;
    Right=(1./(Aphic*beta)).*(Aphic).*(1./(Aphi*beta))';
    
    % Compute the first derivative (with respect to x) of interpolation matrix and its rational approximation at test points
    dxRBFe=dxrbf(rle,c,Yle(:,1)-Yl(:,1)');
    dxAphic=dxRBFe;
    dxRight=(1./(Aphic*beta)).^2.*((dxAphic).*(Aphic*beta)-(Aphic).*(dxAphic*beta)).*(1./(Aphi*beta))';
    
    % Compute the first derivative (with respect to y) of interpolation matrix and its rational approximation at test points
    dyRBFe=dyrbf(rle,c,Yle(:,2)-Yl(:,2)');
    dyAphic=dyRBFe;
    dyRight=(1./(Aphic*beta)).^2.*((dyAphic).*(Aphic*beta)-(Aphic).*(dyAphic*beta)).*(1./(Aphi*beta))';
    
    % Compute the Laplacian operator
    d2RBFe=d2rbf(rle,c);
    d2Aphic=d2RBFe;
    
    d2Right=(1./(Aphic*beta)).^(3).*(((d2Aphic).*(Aphic*beta)-(d2Aphic*beta).*(Aphic)).*(Aphic*beta)-...
        2*(dxAphic*beta).*((dxAphic).*(Aphic*beta)-(Aphic).*(dxAphic*beta))-2*(dyAphic*beta).*((dyAphic).*(Aphic*beta)-(Aphic).*(dyAphic*beta)))...
        .*(1./(Aphi*beta))';
    
    switch (Order)
        case ('0')
            % Local fit
            localfit=Right/AphiR;
        case ('1x')
            % Local fit
            localfit=dxRight/AphiR;
            
        case ('1y')
            % Local fit
            localfit=dyRight/AphiR;
        case ('xx+yy')
            % Local fit
            localfit=d2Right/AphiR;
            
        otherwise
            error('Error')
    end  
    % Assemble 
    A(indle,indl)=A(indle,indl)+localfit.*SEM(indle,l);
end
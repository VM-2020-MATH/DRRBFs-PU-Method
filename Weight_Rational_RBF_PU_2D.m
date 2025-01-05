function A=Weight_Rational_RBF_PU_2D(X,Xe,a,b,n,Order,Tfunc)
% Inputs: X is a set of trial points; 
%         Xe is a set of test points;
%         a and b are the beginning and end of the interval (Domain is the square [a,b]^2)
%         n=(b-a)/h for choosing the number of patch (Npu=(n/4)^2)
%         Order: '1' for the approximation; 
%                '1x' for the first derivative in direction x;
%                '1y' for the first derivative in direction y;
%         Tfunc: is the test function, i.e., 'R1' (first example) or 'R2' (second example)
%         given in the paper;

% Output: A is the differentiation matrix;
         
% Test function
switch (Tfunc)
    case ('R1')
    u=@(xx,yy) atan(125*(sqrt((xx-1.5).^2+(yy-0.25).^2)-0.92));
    case ('R2')
     u=@(xx,yy) tan(9*(yy-xx)+1)./(tan(9)+1);
    otherwise
        error ('No example is found, but the user can add a new example...')
end
% RBF (Matern kernel) and its first derivatives
rbf=@(r,c)       exp(-c*r).*(15+15*(c*r)+6*(c*r).^2+(c*r).^3);
dxrbf=@(r,c,x)  -c^2.*x.*exp(-c*r).*((c*r).^2+3*(c*r)+3);
dyrbf=@(r,c,y) -c^2.*y.*exp(-c*r).*((c*r).^2+3*(c*r)+3);
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
delta=1/npu;
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
    
    % Compute the interpolation matrix (using the MDI)
    mu=0.00000001;
    Aphi=RBFI+mu*eye(Nl); 
    f=u(Yl(:,1),Yl(:,2));
    D=diag(f);

    % Cholesky factorization
    L=chol(Aphi,'lower');
    Ainv=L'\(L\eye(Nl));


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
            
        otherwise
            error('Error')
    end  
    % Assemble 
    A(indle,indl)=A(indle,indl)+localfit.*SEM(indle,l);
end

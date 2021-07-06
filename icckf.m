function [xk1k1,Pk1k1,uk1,Uk1,ak1,bk1]=icckf(xkk,Pkk,y,Q,phi,uk,Uk,ak,bk,N)

nx=size(xkk,1);
nz=size(y,1);
nPts0=2*nx;         %%%%%%%%%%%1st Cubature point 
nPts1=2*(nx+nx);    %%%%%%%%%%%2nd Cubature point

%%%%%%%%%%Xk1k%%%%%%%%%%%%%%%%
Xkk=CR(xkk,Pkk); 

Xk1k=ckf_ProssEq(Xkk);                  %%%%Calculate Cubature point 

xk1k=sum(Xk1k,2)/nPts0;

Pk1k=Xk1k*Xk1k'/nPts0-xk1k*xk1k'+Q;

Pxxkk1k=Xkk*Xk1k'/nPts0-xkk*xk1k';

xak1k=[xk1k;xkk];

Pak1k=[Pk1k Pxxkk1k';Pxxkk1k Pkk];

Xak1k=CR(xak1k,Pak1k);

Yk1k=ckf_Mst_new(Xak1k,phi);

yk1k=sum(Yk1k,2)/nPts1;

Sk=Yk1k*Yk1k'/nPts1-yk1k*yk1k';

Pxyk1k=Xak1k(1:nx,:)*Yk1k'/nPts1-xk1k*yk1k';

Pxykk1k=Xak1k(nx+1:end,:)*Yk1k'/nPts1-xkk*yk1k';

%%%%%%%%Initialization%%%%%%%%%%
uk1=uk+1;
Uk1=Uk;
ak1=ak+0.5;
bk1=bk;
alfa=1;
beta=1;

%%%%%%%%Calculate initial expectation%%%%%%%%
E_lamda=alfa/beta;
E_inv_R=(uk1-nz-1)*inv(Uk1);
E_v=ak1/bk1;

for i=1:N
    
    %%%%%%%%%%R
    R=inv(E_inv_R)/E_lamda;

    %%%%%%%%%%%
    Pyyk1k=Sk+R;
    
    %%%%%%%%%Smoother%%%%%%
    Kks=Pxykk1k*inv(Pyyk1k);

    xkk1=xkk+Kks*(y-yk1k);

    Pkk1=Pkk-Kks*Pyyk1k*Kks';
    
    %%%%%%%Filter%%%%%%%%%%%%
    Kk1=Pxyk1k*inv(Pyyk1k);     

    xk1k1=xk1k+Kk1*(y-yk1k);

    Pk1k1=Pk1k-Kk1*Pyyk1k*Kk1';

    %%%%%%%%Pxxkk1k1%%%%%%%%%
    Pxxkk1k1=Pxxkk1k-Kks*Pyyk1k*Kk1';   
    
    %%%%%%%%Extend state estimation
    xak1k1=[xk1k1;xkk1];
    Pak1k1=[Pk1k1 Pxxkk1k1';Pxxkk1k1 Pkk1];
    
    %%%%%%%%%%%%%%
    Xak1k1=CR(xak1k1,Pak1k1);
    Dk=(repmat(y,1,nPts1)-ckf_Mst_new(Xak1k1,phi))*(repmat(y,1,nPts1)-ckf_Mst_new(Xak1k1,phi))'/nPts1;
    
    %%%%%%%%%Update parameter
    alfa=(nz+E_v)/2;
    beta=(E_v+trace(Dk*E_inv_R))/2;
    
    %%%%%%%%%%Expectation
    E_lamda=alfa/beta;
    E_log_lamda=psi(alfa)-log(beta);
    E_inv_R=(uk1-nz-1)*inv(Uk1);
    E_v=ak1/bk1;
    
    %%%%%%%%%%Update parameter
    Uk1=Uk+E_lamda*Dk;
    bk1=bk-0.5-0.5*E_log_lamda+0.5*E_lamda;

end
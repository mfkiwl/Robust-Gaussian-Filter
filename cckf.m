function [xk1k1,Pk1k1,xkk1,Pkk1,Ak] = cckf(xkk,Pkk,y,Q,R,phi)
nx=size(xkk,1);
nPts0=2*nx;         %%%%%%%%%%%1st Cubature point 
nPts1=2*(nx+nx);    %%%%%%%%%%%2nd Cubature point 

%%%%%%%%%%Xk1k%%%%%%%%%%%%%%%%
Xkk=CR(xkk,Pkk); 

Xk1k=ckf_ProssEq(Xkk);                  %%%%Calculate Cubature point 

xk1k=sum(Xk1k,2)/nPts0;

Pk1k=Xk1k*Xk1k'/nPts0-xk1k*xk1k'+Q;

%%%%%%%%%%%%Pxxkk1k
Pxxkk1k=Xkk*Xk1k'/nPts0-xkk*xk1k';

%%%%%%%%%%%%Pxykk1k
xak1k=[xk1k;xkk];

Pak1k=[Pk1k Pxxkk1k';Pxxkk1k Pkk];

Xak1k=CR(xak1k,Pak1k);

Yk1k=ckf_Mst_new(Xak1k,phi);

yk1k=sum(Yk1k,2)/nPts1;

Pyyk1k=Yk1k*Yk1k'/nPts1-yk1k*yk1k'+R;

Pxyk1k=Xak1k(1:nx,:)*Yk1k'/nPts1-xk1k*yk1k';

Pxykk1k=Xak1k(nx+1:end,:)*Yk1k'/nPts1-xkk*yk1k';

%%%%%%%%%Smoother%%%%%%
Kks=Pxykk1k*inv(Pyyk1k);

xkk1=xkk+Kks*(y-yk1k);

Pkk1=Pkk-Kks*Pyyk1k*Kks';

%%%%%%%Filter%%%%%%%%%%%%
Kk1=Pxyk1k*inv(Pyyk1k);     

xk1k1=xk1k+Kk1*(y-yk1k);

Pk1k1=Pk1k-Kk1*Pyyk1k*Kk1';

%%%%%%%Smoother gain%%%%%%%%
Pxxkk1k1=Pxxkk1k-Kks*Pyyk1k*Kk1';   

Ak=Pxxkk1k1*inv(Pk1k1);

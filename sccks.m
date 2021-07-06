function [xkNA,PkNA]=sccks(xs,Ps,yA,Q,phi,ts,R)

%%%%Initialization
nx=size(xs,1);

%%%%Initial state
xkk=xs;
Pkk=Ps;

%%%%Store state
xkkA=xkk;
PkkA=Pkk;
xk_1kA=[];
Pk_1kA=[];
KsA=[];

for t=1:ts
    
    y=yA(:,t);
        
    %%%%%%Filter process
    [xkk,Pkk,xk_1k,Pk_1k,Ks] = cckf(xkk,Pkk,y,Q,R,phi);

    %%%%Store state
    xkkA=[xkkA xkk];
    PkkA=[PkkA Pkk];
    xk_1kA=[xk_1kA xk_1k];
    Pk_1kA=[Pk_1kA Pk_1k];
    KsA=[KsA Ks];

end

%%%%Initial state
xkN=xkk;
PkN=Pkk;
    
%%%%Store state
xkNA=xkN;
PkNA=PkN;

for t=(ts-1):-1:0
        
    xkk=xkkA(:,t+2);
    Pkk=PkkA(:,(t+1)*nx+1:(t+2)*nx);
    xkk1=xk_1kA(:,t+1);
    Pkk1=Pk_1kA(:,t*nx+1:(t+1)*nx);
    Ak=KsA(:,t*nx+1:(t+1)*nx);

    %%%%%%Smoother process
    [xkN,PkN]=ccks(xkN,PkN,xkk,Pkk,xkk1,Pkk1,Ak);

    %%%%Store state
    xkNA=[xkN xkNA];
    PkNA=[PkN PkNA];

end

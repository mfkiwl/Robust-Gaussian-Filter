function [xkNA,PkNA]=mcks(xi,Pi,zA,ts,Q,R0,v,N)

%%%%Initialization
nx=size(xi,1);
nz=size(zA,1);
E_lamda=ones(1,ts);

for i=1:N

    %%%%Initial state
    xkk=xi;
    Pkk=Pi;
    %%%%Store state
    xkkA=xkk;
    PkkA=Pkk;
    xkk_1A=[];
    Pkk_1A=[];
    Pk_1kk_1A=[]; 
    
    for t=1:ts
        
        %%%%%R
        R=R0/E_lamda(t);
        
        %%%%%%Filter process
        [xkk,Pkk,xkk_1,Pkk_1,Pk_1kk_1]=ckf(xkk,Pkk,zA(:,t),Q,R);   
        
        %%%%Store state
        xkkA=[xkkA xkk];
        PkkA=[PkkA Pkk];
        xkk_1A=[xkk_1A xkk_1];
        Pkk_1A=[Pkk_1A Pkk_1];
        Pk_1kk_1A=[Pk_1kk_1A Pk_1kk_1];

    end
    
    %%%%Initial state
    xkN=xkk;
    PkN=Pkk;
    
    %%%%Store state
    xkNA=xkN;
    PkNA=PkN;

    for t=(ts-1):-1:0
        
        xkk=xkkA(:,t+1);
        Pkk=PkkA(:,t*nx+1:(t+1)*nx);
        xkk_1=xkk_1A(:,t+1);
        Pkk_1=Pkk_1A(:,t*nx+1:(t+1)*nx);
        Pk_1kk_1=Pk_1kk_1A(:,t*nx+1:(t+1)*nx);
        
        %%%%%%Smoother process
        [xkN,PkN,Ks]=cks(xkN,PkN,xkk,Pkk,xkk_1,Pkk_1,Pk_1kk_1);
        
        xkNA=[xkN xkNA];
        PkNA=[PkN PkNA];

    end
    
    for t=1:ts
        
        %%%%%%%%%%%Obtain parameter
        xkN=xkNA(:,t+1);
        PkN=PkNA(:,t*nx+1:(t+1)*nx);
        z=zA(:,t);
        
        %%%%%%%%%%%Auxiliary variable(See Reference)
        XkN=CR(xkN,PkN);
        F_R=(repmat(z,1,2*nx)-ckf_Mst(XkN))*(repmat(z,1,2*nx)-ckf_Mst(XkN))'/(2*nx);
        
        %%%%%%%%%%%lamda
        gama=trace(F_R*inv(R0));
        E_lamda(t)=(v+nz)/(v+gama);

    end

end

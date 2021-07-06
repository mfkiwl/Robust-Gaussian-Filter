function [xkNA,PkNA]=iccks(xi,Pi,yA,Q,phi,ts,a0,b0,u0,U0,N)

%%%%Initialization
nx=size(xi,1);
nz=size(yA,1);

%%%%Initial expectation
%%%%%Measurement model
E_lamda=ones(1,ts);
E_v=a0/b0;
E_IR=(u0-nz-1)*inv(U0);

ui=u0+ts;
ai=a0+0.5*ts;

for i=1:N

    %%%%Initial state
    xkk=xi;
    Pkk=Pi;
    %%%%Store state
    xkkA=xkk;
    PkkA=Pkk;
    xk_1kA=[];
    Pk_1kA=[];
    KsA=[];
    
    for t=1:ts

        %%%%%R
        R=inv(E_IR)/E_lamda(t);
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

        xkNA=[xkN xkNA];
        PkNA=[PkN PkNA];

    end
    
    %%%%Update parameter
    Ui=U0;
    bi=b0-0.5*ts;

    for t=1:ts
        
        %%%%%%%%%%%Obtain parameter
        Gk_1=KsA(:,(t-1)*nx+1:t*nx);
        xk_1N=xkNA(:,t);
        Pk_1N=PkNA(:,(t-1)*nx+1:t*nx);
        %%%%%%%%%%%
        xkN=xkNA(:,t+1);
        PkN=PkNA(:,t*nx+1:(t+1)*nx);
        y=yA(:,t);

        %%%%%%%%%%%Auxiliary parameters(See Reference)
        Pk_1kN=Gk_1*PkN;
        
        %%%%%%%%%%%Auxiliary variable
        [F_R]=FR(xk_1N,xkN,Pk_1N,Pk_1kN,PkN,y,phi);
        
        %%%%%%%%%%%Update parameter of lamda
        alfa_kk=0.5*(nz+E_v);
        beta_kk=0.5*(E_v+trace(F_R*E_IR));
        
        %%%%%%%%%%%E_lamda
        E_lamda(t)=alfa_kk/beta_kk;
        
        %%%%%%%%%%Ui,bi
        Ui=Ui+E_lamda(t)*F_R;
        bi=bi+0.5*E_lamda(t)-0.5*(psi(alfa_kk)-log(beta_kk));
        
    end

    %%%%%%%%%%%%Calculate expectation
    E_v=ai/bi;
    E_IR=(ui-nz-1)*inv(Ui);

end
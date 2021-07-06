%%%%%Preparation%%%%%%%%%%%%%%%%
clear all;
close all;
clc;
randn('state',sum(100*clock));     
format long;

%%%%%Model parameters%%%%%%%%%%%%%%%%
raddeg=pi/180;                     
nx=5;                              
nz=2;                              
nxp=50;                           %%%Monte Carlo simulation number
ts=100;                            
T=1;
M=[T^3/3 T^2/2;T^2/2 T];
%%%%Observation noise%%%%%%%%%%%%%%
q1=0.1;
q2=1.75e-4;
Cr=10;
Co=sqrt(10)*1e-3;
Q=[q1*M zeros(2,2) zeros(2,1);zeros(2,2) q1*M zeros(2,1);zeros(1,2) zeros(1,2) q2*T]; 
R0=diag([Cr^2 Co^2]);

%%%%related parameters
phi=[0.5 0;0 0.5];

%%%%Normal probability
p=0.95;

%%%%Iteration number for variation
N=5;

for expt = 1:nxp
    
    fprintf('MC Run in Process = %d\n',expt); 
    
    %%%%%Initial state%%%%%%%%%%%
    x=[1000;300;1000;0;-3*raddeg];        %%%Real initial value
    P=diag([100 10 100 10 1e-4]);         %%%Initial P
    Skk=utchol(P);                        
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%Generate initial measurement
    test=rand;
    if test<=p
        R=R0;
    else
        R=100*R0;
    end
    %%%%Calculate root mean square matrix
    SR=utchol(R);
    v=SR*randn(nz,1);
    z=MstEq(x)+v;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%Initial value for improved CKS
    xi=x+Skk*randn(nx,1);                 %%%Initial state
    Pi=P;
    
    %%%%Initial state of MCKS
    xm=xi;
    Pm=Pi;
    
    %%%%Initial state of CKS
    xs=xi;
    Ps=Pi;
    
    %%%%Initial state of CKF
    xif=xi;
    Pif=Pi;
    aif=5;
    bif=1;
    uif=nz+2;
    Uif=R0;
    
    %%%%Initial state of MCKF
    xmf=xi;
    Pmf=Pi;
    
    %%%%Initial state of colored CKF
    xsf=xi;
    Psf=Pi;

    %%%%Store data
    xA=x;
    zA=z;
    yA=[];
    xifA=xif;
    xmfA=xmf;
    xsfA=xsf;

    %%%%Simulation real state and measurement
    for t=1:ts
        
        test=rand;

        if test<=p
            R=R0;
        else
            R=100*R0;
        end

        %%%%Calculate mean root square matrix
        SQ=utchol(Q);    
        SR=utchol(R);
        
        %%%%Simulation real state and measurement
        x=ProssEq(x)+SQ*randn(nx,1);
        %%%%Generate colored heavy-tailed noise
        v=phi*v+SR*randn(nz,1);
        %%%%Generate measurement
        z=MstEq(x)+v;
        
        %%%%Generate white measurement
        y=z-phi*zA(:,t);
        
        %%%%Filtering process
        [xif,Pif,uif,Uif,aif,bif]=icckf(xif,Pif,y,Q,phi,uif,Uif,aif,bif,N);
        
        [xmf,Pmf]=orckf(xmf,Pmf,z,Q,R0,5,N); 
        
        [xsf,Psf,xsfkk1,Psfkk1,Ask]=cckf(xsf,Psf,y,Q,R0,phi);

        %%%%Store state and measurement
        xA=[xA x];  
        zA=[zA z];
        yA=[yA y];
        
        xifA=[xifA xif];  
        xmfA=[xmfA xmf];  
        xsfA=[xsfA xsf];  
        
    end

    %%%%Initial parameters for improved CKS
    a0=5;
    b0=1;
    u0=nz+2;
    U0=R0;
    %%%%Obtain measurements
    zA=zA(:,2:end);
    
    %%%%Smoother process
    [ixsBA,iPsBA]=iccks(xi,Pi,yA,Q,phi,ts,a0,b0,u0,U0,N);
    
    [mxsBA,mPsBA]=mcks(xm,Pm,zA,ts,Q,R0,5,N);
    
    [sxsBA,sPsBA]=sccks(xs,Ps,yA,Q,phi,ts,R0);

    %%%%Calculate MSE
    %%%%%%Filter
    mse_ickf_1(expt,:)=(xA(1,:)-xifA(1,:)).^2+(xA(3,:)-xifA(3,:)).^2;
    mse_ickf_2(expt,:)=(xA(2,:)-xifA(2,:)).^2+(xA(4,:)-xifA(4,:)).^2;
    mse_ickf_3(expt,:)=(xA(5,:)-xifA(5,:)).^2;
    
    mse_mckf_1(expt,:)=(xA(1,:)-xmfA(1,:)).^2+(xA(3,:)-xmfA(3,:)).^2;
    mse_mckf_2(expt,:)=(xA(2,:)-xmfA(2,:)).^2+(xA(4,:)-xmfA(4,:)).^2;
    mse_mckf_3(expt,:)=(xA(5,:)-xmfA(5,:)).^2;
    
    mse_sckf_1(expt,:)=(xA(1,:)-xsfA(1,:)).^2+(xA(3,:)-xsfA(3,:)).^2;
    mse_sckf_2(expt,:)=(xA(2,:)-xsfA(2,:)).^2+(xA(4,:)-xsfA(4,:)).^2;
    mse_sckf_3(expt,:)=(xA(5,:)-xsfA(5,:)).^2;
    
    %%%%%Smoother
    mse_icks_1(expt,:)=(xA(1,:)-ixsBA(1,:)).^2+(xA(3,:)-ixsBA(3,:)).^2;
    mse_icks_2(expt,:)=(xA(2,:)-ixsBA(2,:)).^2+(xA(4,:)-ixsBA(4,:)).^2;
    mse_icks_3(expt,:)=(xA(5,:)-ixsBA(5,:)).^2;
    
    mse_mcks_1(expt,:)=(xA(1,:)-mxsBA(1,:)).^2+(xA(3,:)-mxsBA(3,:)).^2;
    mse_mcks_2(expt,:)=(xA(2,:)-mxsBA(2,:)).^2+(xA(4,:)-mxsBA(4,:)).^2;
    mse_mcks_3(expt,:)=(xA(5,:)-mxsBA(5,:)).^2;
    
    mse_scks_1(expt,:)=(xA(1,:)-sxsBA(1,:)).^2+(xA(3,:)-sxsBA(3,:)).^2;
    mse_scks_2(expt,:)=(xA(2,:)-sxsBA(2,:)).^2+(xA(4,:)-sxsBA(4,:)).^2;
    mse_scks_3(expt,:)=(xA(5,:)-sxsBA(5,:)).^2;

end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%Filter
rmse_ickf_1=sqrt(mean(mse_ickf_1,1));
rmse_ickf_2=sqrt(mean(mse_ickf_2,1));
rmse_ickf_3=sqrt(mean(mse_ickf_3,1));

rmse_mckf_1=sqrt(mean(mse_mckf_1,1));
rmse_mckf_2=sqrt(mean(mse_mckf_2,1));
rmse_mckf_3=sqrt(mean(mse_mckf_3,1));

rmse_sckf_1=sqrt(mean(mse_sckf_1,1));
rmse_sckf_2=sqrt(mean(mse_sckf_2,1));
rmse_sckf_3=sqrt(mean(mse_sckf_3,1));

%%%%%%%Smoother
rmse_icks_1=sqrt(mean(mse_icks_1,1));
rmse_icks_2=sqrt(mean(mse_icks_2,1));
rmse_icks_3=sqrt(mean(mse_icks_3,1));

rmse_mcks_1=sqrt(mean(mse_mcks_1,1));
rmse_mcks_2=sqrt(mean(mse_mcks_2,1));
rmse_mcks_3=sqrt(mean(mse_mcks_3,1));

rmse_scks_1=sqrt(mean(mse_scks_1,1));
rmse_scks_2=sqrt(mean(mse_scks_2,1));
rmse_scks_3=sqrt(mean(mse_scks_3,1));

%%%%%%%%%%%Plot
figure;
j = 0:ts;
plot(j,rmse_sckf_1(1,:),'--g',j,rmse_mckf_1(1,:),'--b',j,rmse_ickf_1(1,:),'--r',j,rmse_scks_1(1,:),'-g',j,rmse_mcks_1(1,:),'-b',j,rmse_icks_1(1,:),'-r','linewidth',2.5);
xlabel('Time (s)');
ylabel('Position RMSE (m)');
legend('Colored CKF','Robust CKF','Improved robust CKF','Colored CKS','Robust CKS','Improved robust CKS');

figure;
j = 0:ts;
plot(j,rmse_sckf_2(1,:),'--g',j,rmse_mckf_2(1,:),'--b',j,rmse_ickf_2(1,:),'--r',j,rmse_scks_2(1,:),'-g',j,rmse_mcks_2(1,:),'-b',j,rmse_icks_2(1,:),'-r','linewidth',2.5);
xlabel('Time (s)');
ylabel('Velocity RMSE (m/s)');
legend('Colored CKF','Robust CKF','Improved robust CKF','Colored CKS','Robust CKS','Improved robust CKS');

figure;
j = 0:ts;
plot(j,rmse_sckf_3(1,:)./raddeg,'--g',j,rmse_mckf_3(1,:)./raddeg,'--b',j,rmse_ickf_3(1,:)./raddeg,'--r',j,rmse_scks_3(1,:)./raddeg,'-g',j,rmse_mcks_3(1,:)./raddeg,'-b',j,rmse_icks_3(1,:)./raddeg,'-r','linewidth',2.5);
xlabel('Time (s)');
ylabel('Turn Rate RMSE (Deg/s)');
legend('Colored CKF','Robust CKF','Improved robust CKF','Colored CKS','Robust CKS','Improved robust CKS');


%%%%%%%%%%%%%%%%%%
armse_ickf_1=sqrt(mean(mean(mse_ickf_1,1)))
armse_ickf_2=sqrt(mean(mean(mse_ickf_2,1)))
armse_ickf_3=sqrt(mean(mean(mse_ickf_3,1)))./raddeg

armse_icks_1=sqrt(mean(mean(mse_icks_1,1)))
armse_icks_2=sqrt(mean(mean(mse_icks_2,1)))
armse_icks_3=sqrt(mean(mean(mse_icks_3,1)))./raddeg

armse_mckf_1=sqrt(mean(mean(mse_mckf_1,1)))
armse_mckf_2=sqrt(mean(mean(mse_mckf_2,1)))
armse_mckf_3=sqrt(mean(mean(mse_mckf_3,1)))./raddeg
 
armse_mcks_1=sqrt(mean(mean(mse_mcks_1,1)))
armse_mcks_2=sqrt(mean(mean(mse_mcks_2,1)))
armse_mcks_3=sqrt(mean(mean(mse_mcks_3,1)))./raddeg

armse_sckf_1=sqrt(mean(mean(mse_sckf_1,1)))
armse_sckf_2=sqrt(mean(mean(mse_sckf_2,1)))
armse_sckf_3=sqrt(mean(mean(mse_sckf_3,1)))./raddeg
 
armse_scks_1=sqrt(mean(mean(mse_scks_1,1)))
armse_scks_2=sqrt(mean(mean(mse_scks_2,1)))
armse_scks_3=sqrt(mean(mean(mse_scks_3,1)))./raddeg


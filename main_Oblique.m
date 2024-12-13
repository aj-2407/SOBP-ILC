clear;
clc;
close all;

%oblique projection parameter setting
nb=10;
dem=nb+1;
Ts = 0.0005;%sampling time
T_end = 9-Ts;
TIMES=11;%itertion index
tt=1;
t = 0:tt*Ts:T_end;
t =t';
N = length(t);
one(1:N,1)=1;

%initialization
theta_true=zeros(dem,1);theta=zeros(dem,TIMES);
err=zeros(N,TIMES);delta=zeros(dem,TIMES); 
TrackErrNorm=zeros(1,TIMES);ParamErrNorm=zeros(dem,TIMES);
D3=zeros(nb+1,3600);




load('r.mat');%REFERENCE
load('dr.mat');%REFERENCE VELOCITY
load('gb_18000.mat');%spatial_impulse matrix
Gb = zeros(N,N);
for i = 1:N
    Gb(i:N,i) = gb(1:N+1-i); % Toeplitz matrix
end  
             

%Matrix based on Bernstein polynomials 
for jj=1:nb+1
    for ii=1:N
        D3(jj,ii)=(ii/N)^(jj-1)*(1-ii/N)^(nb+1-jj);
    end
end
Ds3=D3';
PHI=Ds3;
%QR decomposition
[Q,R] = QR_decompos(PHI);  
A = Q;
%Obtain oblique projection matrices
M=A'*Gb*PHI;


 
%identification period
for j=1:TIMES-1
    %start iteration loop

    if j<=20
    Uff=PHI*theta(:,j);
    simin_r=[t r];%Change the input signal
    simin_Uff=[t Uff];%Change the feedforward signal
    sim('sys_360'); %run the simulink
    err(:,j)=r-y;

    %Project error signal based on the matrices according to (19) and (24)
    delta(:,j)=A'*(err(:,j));
    %Update parameter estimates
    theta(:,j+1)=theta(:,j)-0.7*inv(M)*delta(:,j);
    end

    
    %Obtain ErrNorm
    TrackErrNorm(j)=norm(err(:,j));
    ParamErrNorm(j)=norm(theta(:,j)-theta_true);
    
    
end

%set(0,'defaultfigurecolor','w');
figure(1);
plot(t,err(:,1),'r','LineWidth',2)%1st
legend('1st without ILC');
xlabel('Time (s)');ylabel('Tracking error (deg)')
ylim([-0.2,0.2]);

figure(2);
plot(t,err(:,TIMES-1),'b','LineWidth',2);hold on;%9th
xlabel('Time (s)');ylabel('Tracking error (deg)');
legend('10th SOBP-ILC Case 2');

figure(3);
plot(t,err(:,1),'r','LineWidth',2);hold on;%1st
plot(t,err(:,2),'Color','[0.93,0.69,0.13]','LineWidth',2);%,'linestyle','-.');
plot(t,err(:,TIMES-1),'b','LineWidth',2);%,'linestyle','--');
legend('1st without ILC','2nd SOBP-ILC Case 2','10th SOBP-ILC Case 2');
xlabel('Time (s)');ylabel('Tracking error (deg)');
ylim([-0.2,0.2]);

figure(4);
x = 1:1:TIMES;
X=1:1:TIMES;
y1=spline(x,TrackErrNorm,X);plot(X,y1,'r','LineWidth',2);
xlabel('Iteration');ylabel('||Error||_2');
xlim([1,10]);


figure(5)
plot(t,Uff,'LineWidth',2); hold on;
xlabel('Time (s)');ylabel('Feedforward signal (V)');
ylim([-2,3]);

% figure(6);
% for jj=1:nb+1
%     plot(theta(jj,:),'LineWidth',2); hold on;%,'Color','[0,0.45,0.74]'
% end
% xlim([1,TIMES]);
% xlabel('Iteration');ylabel('Estimates');
 


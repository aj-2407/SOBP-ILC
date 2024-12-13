clear;
clc;
close all;


Ts = 0.0005;
T_end = 9-Ts;
TIMES=11;
dem=3;


tt=1;
t = 0:tt*Ts:T_end;
t =t';
N = length(t);
one(1:N,1)=1;
theta_true=zeros(dem,1);
theta=zeros(dem,TIMES);
err=zeros(N,TIMES);TrackErrNorm=zeros(1,TIMES);ParamErrNorm=zeros(dem,TIMES);
delta=zeros(dem,TIMES); 

load('r.mat');%REFERENCE
load('dr.mat');%REFERENCE VELOCITY
load('gb_18000.mat');%spatial_impulse matrix

Gb = zeros(N,N);
for i = 1:N
    Gb(i:N,i) = gb(1:N+1-i); % Toeplitz¾ØÕó
end  


% %temporal --> spatial
% 
% for i1=1:N
%     for i=2:N*5
%         if fix(10*r(i))==i1-1%fix(r_w(i))==i1%mod(y_spatial(i),i1)==i1
%             r_w2(i1,:)=r(i-1);%r_w2(i1+1,:)=r_w(i-1);
%             dr_w2(i1,:)=dr(i-1);ddr_w2(i1,:)=ddr(i-1);
%             break;
%         end
%     end
% end





%PHI feedforward signal using COP-ILC basis functions
%PHI=[dr(1:N),-sin(2*pi/360*r(1:N)),-cos(2*pi/360*r(1:N)),-sin(pi/12*r(1:N)),-cos(pi/12*r(1:N)),-sin(pi/5*r(1:N)),-cos(pi/5*r(1:N))];
PHI=[-sin(2*pi/360*r(1:N)),-cos(2*pi/360*r(1:N)),dr];



% projection matrix
[Q,R] = QR_decompos(PHI);  %QR decomposition
A = Q;
M=A'*Gb*PHI;

 
%START
for j=1:TIMES-1

    if j<=20
    Uff=PHI*theta(:,j);
    simin_r=[t r];
    simin_Uff=[t Uff];
    sim('sys_360');
    err(:,j)=r-y;

    
    delta(:,j)=A'*(err(:,j));
    theta(:,j+1)=theta(:,j)-0.7*inv(M)*delta(:,j);

    end

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
plot(t,err(:,TIMES-1),'b','LineWidth',2);%9th
xlabel('Time (s)');ylabel('Tracking error (deg)');
legend('10th COP-ILC Case 1');

figure(3);
plot(t,err(:,1),'r','LineWidth',2);hold on;%1st
plot(t,err(:,2),'Color','[0.93,0.69,0.13]','LineWidth',2);%,'linestyle','-.');
plot(t,err(:,TIMES-1),'b','LineWidth',2);%,'linestyle','--');
legend('1st without ILC','2nd COP-ILC Case 1','10th COP-ILC Case 1');
xlabel('Time (s)');ylabel('Tracking error (deg)');
ylim([-0.2,0.2]);

figure(4);
% x = 1:1:20;
% X=1:0.001:20;
% y1=spline(x,TrackErrNorm,X);
% plot(X,y1,'r','LineWidth',0.7);hold on;
plot(TrackErrNorm,'r','LineWidth',2); hold on;
xlim([1,10]);legend('COP-ILC Case 1');
xlabel('Iteration');ylabel('||Error||_2');


figure(6);
for jj=1:dem
    plot(theta(jj,:),'LineWidth',2); hold on;%,'Color','[0,0.45,0.74]'
end
xlim([1,10]);%legend('a','b','c');
xlabel('Iteration');ylabel('Estimates');
clc
clear
close all

global R d dt

% State and Initialize
dt=0.001;
tf=1;
t=[0:dt:tf]';
m=length(t);

N_filter_bank=100; %No. of parallel filters

h=[1 0 0;0 1 0];
r=0.0001*eye(2);

R=0.05; %wheel radius
d=.1; %half of axle length
u=[50 40]';

x0=[0 0 0]';
n=length(x0);
ym=zeros(2,m);

p0=0.001*eye(n); %we know our initial condition of the states with high confidence
p=p0; %changes several times

pe=kron(ones(1,N_filter_bank),p0); %state covariances
xe=kron(ones(1,N_filter_bank),x0); %state estimates

p_save=zeros(3,m);
p_save(:,1)=diag(p0);

% Process Noise (note: there is coupling but is ignored)
q_truth=.005;
qmin=.001;
qmax=.020;
q=q_truth*eye(3);


% Range of values we guess for process noise covariance
q_filter_bank=linspace(qmin,qmax,N_filter_bank)';

%initial uniform weight of filters
w=ones(N_filter_bank,1)/N_filter_bank;

%we update our estimate of q, m times. m=size(measurement samples)
est_q=zeros(m,1);
est_q(1)=q_filter_bank'*w; %initial estimate of q

err_est_q=1000*ones(N_filter_bank,m);
cov_est_q=1000*ones(m,1);

err_est_q(:,1)=q_filter_bank-est_q(1);
cov_est_q(1)=err_est_q(:,1)'*(err_est_q(:,1).*w);

[t,ym]=generate_measurement(t,x0,u,h,r,q_truth);

x_est=zeros(3,length(t));
x_est(:,1)=x0;

% est_cov_x=zeros(3,3*N_filter_bank);
sig2_est_x=zeros(3,m);


for i=1:m-1
    
    ly_given_x=zeros(N_filter_bank,1);
    
    for j=1:N_filter_bank
    
        qj=q_filter_bank(j);
        xej=xe(:,j);
        p=pe(:,(j-1)*n+1:j*n);
        
        % Kalman Update
        r_cov=h*p*h'+r;
        gain=p*h'*inv(r_cov);
        p=(eye(3)-gain*h)*p;
        xej=xej+gain*(ym(:,i)-h*xej);
        
        %Propagation
        f1=dt*diff_drive(xej,u);f2=dt*diff_drive(xej+0.5*f1,u);
        f3=dt*diff_drive(xej+0.5*f2,u);f4=dt*diff_drive(xej+f3,u);
        [Fc,bc]=linearize1(xej,u);
        phi=c2d(Fc,bc,dt);
        
        
        p=phi*p*phi'+qj*dt;
        pe(:,(j-1)*n+1:j*n)=p;
        
        tmp_err=ym(:,i)-h*xej;
        ly_given_x(j)=det(2*pi*r_cov)^(-.5) *exp(-.5*tmp_err'*inv(r_cov)*tmp_err);

        xej=xej+1/6*(f1+2*f2+2*f3+f4);
        xe(:,j)=xej; 
    end
    
    w=w.*ly_given_x;
    w=w/sum(w);
    
    est_q(i+1)=q_filter_bank'*w;
    x_est(:,i+1)=sum(xe'.*w)';
    
    %xe,3by100 pe,3by300 x_est(:,i+1)3by1
    tmp=zeros(n,n);
    for j=1:N_filter_bank
       tmp=tmp+w(j)*(   (xe(:,j)-x_est(:,i+1))*(xe(:,j)-x_est(:,i+1))'    +pe(:,(j-1)*n+1:j*n));
    end

    sig2_est_x(:,i)=diag(tmp);
   
    err_est_q(:,i+1)=q_filter_bank-est_q(i+1);
    cov_est_q(i+1)=err_est_q(:,i+1)'*(err_est_q(:,i+1).*w);

end

% 3-Sigma Outlier
sig3=cov_est_q.^(0.5)*3;
figure;subplot(211);plot(t,sig3,'k',t,-sig3,'k',t,(est_q-q_truth),'b');grid on;
subplot(212);plot(q_filter_bank,w);grid on;

figure;
subplot(122);plot(x_est(1,:),x_est(2,:),'k');grid on;
subplot(121);plot(ym(1,:),ym(2,:),'*');grid on;


[t,x_det]=ode45(@(t,x) f(t,x,u,0,0),t,x0);

% use the estimated q in an EKF to check the result.
sig3_est_x=zeros(3,m);
p0=.001*eye(3);

pe=kron(ones(1,m),p0); %state covariances
xe=zeros(n,m);
xe(:,1)=x0;

qj=est_q(end);

for i=1:m-1
        
        xej=xe(:,i);
        p=pe(:,(i-1)*n+1:i*n);
        
        % Kalman Update
        r_cov=h*p*h'+r;
        gain=p*h'*inv(r_cov);
        p=(eye(3)-gain*h)*p;
        xej=xej+gain*(ym(:,i)-h*xej);
        
        %Propagation
        f1=dt*diff_drive(xej,u);f2=dt*diff_drive(xej+0.5*f1,u);
        f3=dt*diff_drive(xej+0.5*f2,u);f4=dt*diff_drive(xej+f3,u);
        [Fc,bc]=linearize1(xej,u);
        phi=c2d(Fc,bc,dt);

        p=phi*p*phi'+qj*dt;
        pe(:,i*n+1:(i+1)*n)=p;
        sig3_est_x(:,i+1)=diag(p);
        xej=xej+1/6*(f1+2*f2+2*f3+f4);
        xe(:,i+1)=xej; 
end

figure;
subplot(311);plot(t,xe(1,:)-x_det(:,1)',t,3*sqrt(sig3_est_x(1,:)),...
    'k',t,-3*sqrt(sig3_est_x(1,:)));set(gca, 'fontsize',20);
subplot(312);plot(t,xe(2,:)-x_det(:,2)',t,3*sqrt(sig3_est_x(2,:)),...
    'k',t,-3*sqrt(sig3_est_x(2,:)));set(gca, 'fontsize',20);
subplot(313);plot(t,xe(3,:)-x_det(:,3)',t,3*sqrt(sig3_est_x(3,:)),...
    'k',t,-3*sqrt(sig3_est_x(3,:)));set(gca, 'fontsize',20);

keyboard

function xdot=diff_drive(x,u)
    global R d 
    
    phi=x(3);
    mat=[0.5*R*cos(phi) 0.5*R*cos(phi);0.5*R*sin(phi) 0.5*R*sin(phi);-.5*R/d +.5*R/d];
    xdot=mat*u;
end

function [t,ym]=generate_measurement(tspan,x_0,u,h,r,q_truth)
    m=length(tspan);
    n=length(x_0); %system order
    
    [t,x]=ode45(@(t,x) f(t,x,u,1,q_truth),tspan,x_0); %running sys with process noise
    ym=h*x';
    ym=ym+sqrt(r(1,1))*randn(size(ym)); %obtain measurements by applying measurement noise
end

function xdot=f(t,xx,u,coeff,q)
    
    global R d dt
    %coeff is either 0 or 1. the purpose is sometimes we want to run the
    %system without noise =>> coeff=0
    
    x=xx(1);
    y=xx(2);
    phi=xx(3);

    n=length(xx);

    w=sqrt(dt)*q*randn(n,1);
    
    mat=[0.5*R*cos(phi) 0.5*R*cos(phi);
        0.5*R*sin(phi) 0.5*R*sin(phi);
        -.5*R/d .5*R/d];
    tmp=mat*u;
    xdot(1,1)=tmp(1)+coeff*w(1);
    xdot(2,1)=tmp(2)+coeff*w(2);
    xdot(3,1)=tmp(3)+coeff*w(3);
    
end

function [F,b]=linearize1(x_0,u)
    
    global R d 
    phi=x_0(3);
    uL=u(1);
    uR=u(2);
    
    F=[0, 0, - (R*uL*sin(phi))/2 - (R*uR*sin(phi))/2;
       0, 0,   (R*uL*cos(phi))/2 + (R*uR*cos(phi))/2;
       0, 0 ,0];
    b=[(R*cos(phi))/2, (R*cos(phi))/2;
       (R*sin(phi))/2, (R*sin(phi))/2;
       -R/(2*d),        R/(2*d)];
end
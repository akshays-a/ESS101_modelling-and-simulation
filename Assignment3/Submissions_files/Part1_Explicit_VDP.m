clc
clear
close

%given conditions 
u  = 5; 
x0 = [1;0];
tfinal=25;

%% ODE solvers

%using ODE45 default
xOde45_def  = sym('x',[2,1]);
syms tODE45_def real

f1=[xOde45_def(2); u*(1-xOde45_def(1)^2)*xOde45_def(2)-xOde45_def(1)];
matlabFunction(f1, 'file', 'VanDerPol','vars',{tODE45_def,xOde45_def});

options1=odeset();
[tODE45_def,xOde45_def] = ode45(@VanDerPol, [0 tfinal], x0,options1);

%using ODE45 with tighter tolerance
xOde45_tight  = sym('x',[2,1]);
syms tODE45_tight real

f2=[xOde45_tight(2); u*(1-xOde45_tight(1)^2)*xOde45_tight(2)-xOde45_tight(1)];
matlabFunction(f2, 'file', 'VanDerPol','vars',{tODE45_tight,xOde45_tight});

options2=odeset('AbsTol',1e-8,'RelTol',1e-8);
[tODE45_tight,xOde45_tight] = ode45(@VanDerPol, [0 tfinal], x0,options2);


%calculating the time step for ODE45 solvers
dtOde45_def=zeros(length(tODE45_def),1);
for i= 2: length(tODE45_def)
    dtOde45_def(i-1)=tODE45_def(i)-tODE45_def(i-1);
end

dtOde45_tight=zeros(length(tODE45_tight),1);
for m= 2: length(tODE45_tight)
    dtOde45_tight(m-1)=tODE45_tight(m)-tODE45_tight(m-1);
end
%% RK4 method

f_x = @(yRK4)(yRK4);
f_y = @(xRK4,yRK4)(u*(1-xRK4^2)*yRK4-xRK4);

dt= 10^-2;
N = tfinal/dt;
tRK4 = zeros(N+1,1);
xRK4 = zeros(N,1);
yRK4 = zeros(N,1);
xRK4(1)=x0(1);
yRK4(1)=x0(2);

    for j = 2:N+1
        ButcherRK4=butchert(4);
        a21=ButcherRK4(2,2);
        a32=ButcherRK4(3,3);
        a43=ButcherRK4(4,4);
        c=ButcherRK4(:,1)';
        b=ButcherRK4(5,2:end)';
        tRK4(j) = tRK4(j-1) + dt;
        
        K1x = f_x(yRK4(j-1));
        K1y = f_y(xRK4(j-1), yRK4(j-1));
        K2x = f_x(yRK4(j-1)+a21*dt*K1y);
        K2y = f_y(xRK4(j-1)+a21*dt*K1x, yRK4(j-1)+c(2)*dt*K1y);
        K3x = f_x(yRK4(j-1)+a32*dt*K2y);
        K3y = f_y(xRK4(j-1)+a32*dt*K2x, yRK4(j-1)+c(3)*dt*K2y);
        K4x = f_x(yRK4(j-1)+a43*dt*K3y);
        K4y = f_y(xRK4(j-1)+a43*dt*K3x, yRK4(j-1)+c(4)*dt*K3y);
        xRK4(j) = xRK4(j-1) + dt*(b(1)*K1x+ b(2)*K2x+b(3)*K3x+b(4)*K4x);
        yRK4(j) = yRK4(j-1) + dt*(b(1)*K1y+ b(2)*K2y+b(3)*K3y+b(4)*K4y);
    end
x=ones(100,1);
for k=1:100
    dtRK4(k,1)=dt*x(k);
end

%% Plots

figure(1);clf;
title('Van Der Pol')
subplot(3,3,1)
plot(tODE45_def,xOde45_def(:,1),'r');
hold on
plot(tODE45_def,xOde45_def(:,2),'b'); 
hold on
legend('default ode45')
subplot(3,3,2)
plot(tRK4,xRK4,'b')
hold on
plot(tRK4,yRK4,'r')
xlabel('t'); ylabel('x(t)'); 
legend('xRK4')
xlim([0 25])
subplot(3,3,3)
plot(tODE45_tight,xOde45_tight(:,1),'m');
hold on
plot(tODE45_tight,xOde45_tight(:,2),'c'); 
legend('tighter tol ode45')
subplot(3,3,4)
plot(dtOde45_def,'--r')
legend('default dtode45')
subplot(3,3,5)
plot(dtRK4,'--b')
legend('time step RK4')
subplot(3,3,6)
plot(dtOde45_tight,'--m')
legend('tighter tol dtode45')
hold off
%% 

figure(2);clf
plot(tRK4,xRK4,'or',tODE45_tight,xOde45_tight(:,1),'b')
hold on
plot(tRK4,yRK4,'or',tODE45_tight,xOde45_tight(:,2),'b')
legend('RK4','ODE45 tighter tolerance')
xlim([0 25])
clc;clear;
close all

u=5;
f_x = @(y)(y);
f_y = @(x,y)(u*(1-x^2)*y-x);

%Simulation parameters:
x0=1;
y0=0;
tFinal = 25;


dt= 10^-3;
N = tFinal/dt;
t = zeros(N+1,1);
xRK4 = zeros(N,1);
yRK4 = zeros(N,1);
xRK4(1)=x0;
yRK4(1)=y0;

    for j = 2:N+1

        ButcherRK4=butchert(4);
        a21=ButcherRK4(2,2);
        a33=ButcherRK4(3,3);
        a44=ButcherRK4(4,4);
        c=ButcherRK4(:,1)';
        b=ButcherRK4(5,2:end)';
        t(j) = t(j-1) + dt;
        
        K1x = f_x(yRK4(j-1));
        K1y = f_y(xRK4(j-1), yRK4(j-1));
        K2x = f_x(yRK4(j-1)+a21*dt*K1y);
        K2y = f_y(xRK4(j-1)+a21*dt*K1x, yRK4(j-1)+c(2)*dt*K1y);
        K3x = f_x(yRK4(j-1)+a33*dt*K2y);
        K3y = f_y(xRK4(j-1)+a33*dt*K2x, yRK4(j-1)+c(3)*dt*K2y);
        K4x = f_x(yRK4(j-1)+a44*dt*K3y);
        K4y = f_y(xRK4(j-1)+a44*dt*K3x, yRK4(j-1)+c(4)*dt*K3y);
        xRK4(j) = xRK4(j-1) + dt*(b(1)*K1x+ b(2)*K2x+b(3)*K3x+b(4)*K4x);
        yRK4(j) = yRK4(j-1) + dt*(b(1)*K1y+ b(2)*K2y+b(3)*K3y+b(4)*K4y);

    end

plot(t,xRK4     ,'marker','.','markersize',10)
%hold on
%plot(t,yRK4     ,'marker','.','markersize',10)
xlabel('t');
ylabel('x');
legend('xRK4')

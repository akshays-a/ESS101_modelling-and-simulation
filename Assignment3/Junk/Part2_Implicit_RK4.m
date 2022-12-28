clc
clear
close all

%% Defining the function:
syms lambda x t real

f = lambda*x;

matlabFunction(f, 'file', 'Testfn','vars',{lambda,x});
clear x t f
%% Defining r(xNext,x,u): (It is better to define this in a separate file)
K = sym('K',[2,1]);
syms xk dt real
A= [0.25 0.25-(sqrt(3)/6);0.25+(sqrt(3)/6) 1/4];
a11=A(1 ,1);
a12=A(1 ,2);
a21=A(2 ,1);
a22=A(2 ,2);
b= [0.5 0.5];
tFinal=1;

r=[Testfn(lambda, xk+(dt*a11*K(1)+dt*a12*K(2)))-K(1); Testfn(lambda, xk+(dt*a21*K(1)+dt*a22*K(2)))-K(2)];

dr = jacobian(r,K);

matlabFunction(r,dr, 'file', 'rFileIRK4','vars',{lambda,xk,K,dt});
clear K x t dt r dr
%K=[f(lambda, x+(dt*A(1,1)*K1+dt*A(1,2)*K2));f(lambda, x+(dt*A(2,1)*K1+dt*A(2,2)*K2))];
% for s=1:stages
%     for j=1:stages
%         for i=1:stages
%             K(s)=f(lambda, x+(dt*(A(i,j)*K(j))));
%         end
%     end
% end

%% Simulation:
% Simulate using IEuler
tFinal = 10;
dt = 0.1;
x0 = [1;1];
Nsteps = tFinal/dt;
t = 0:dt:tFinal;
xRK4 = 0;
K = x0;
lambda=-2;
xNext=0;

% Loop for the Implicit Euler
for k = 1:Nsteps
    % Newton iteration
    iter = true;
    alpha = 1;
    niter = 0;
    while iter
        [r,dr] = rFileIRK4(lambda,xRK4,K,dt);
     
        K = K - dr\r;
        
        norm(r);
        if norm(r) < 1e-5
            iter = false;
        else
            niter = niter + 1;
        end
    end
    xNext=xk+(dt*b(1 ,2)*K(1,1)+dt*b(2 ,2)*K(2,1));
    xRK4 = xNext;
end

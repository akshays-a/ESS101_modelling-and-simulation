Stages=2;
n=2;
u=5;
tf=25;
dt=0.01;
N=tf/dt;

%Here we compute the Implicit RK4:
ButchTIM=[1/2 - sqrt(3)/6, 1/4, 1/4-sqrt(3)/6;1/2 + sqrt(3)/6, 1/4 + sqrt(3)/6, 1/4;0, 1/2, 1/2];
aRK4IM=ButchTIM(1:2,2:3);  
bRK4IM=ButchTIM(3,2:3);
cRK4IM=ButchTIM(1:2,1);
KIM_int=zeros(n*Stages,n*Stages)';
xRK4IM=ones(N,1);
yRK4IM=ones(N,1);
KIM=sym('KIM',[n,Stages],'real');
r1=sym('r1',[1,Stages*n],'real')';
KIM=reshape(KIM,1,Stages*n);

%Here compute the necessary functions for the Newton method, and fill in r1:
for j = 1:Stages
    r1(n*j-(n-1):n*j-(n-1))=KIM(n*j-(n-1))-ht1(x+dt*(aRK4IM(j,1)*KIM(1)+aRK4IM(j,2)*KIM(3)),y+dt*(aRK4IM(j,1)*KIM(2)+aRK4IM(j,2)*KIM(4)));
end

for j = 1:Stages
     r1(n*j-(n-1)+1:n*j-(n-1)+1)=KIM(n*j-(n-1)+1)-ht2(x+dt*(aRK4IM(j,1)*KIM(1)+aRK4IM(j,2)*KIM(3)),y+dt*(aRK4IM(j,1)*KIM(2)+aRK4IM(j,2)*KIM(4)));
end
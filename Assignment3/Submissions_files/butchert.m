function [Btab] = butchert(order)

Btab=zeros(order+1,order+1);
if order == 1
    Btab(2,2)=1;
end

if order == 2
    Btab(2,1)=0.5;
    Btab(2,2)=0.5;
    Btab(3,3)=1;
end

if order == 4
    Btab(2,1)=0.5;
    Btab(3,1)=0.5;
    Btab(4,1)=1;
    Btab(2,2)=0.5;
    Btab(3,3)=0.5;
    Btab(4,4)=1;
    Btab(5,2)=1/6;
    Btab(5,3)=1/3;
    Btab(5,4)=1/3;
    Btab(5,5)=1/6;
end

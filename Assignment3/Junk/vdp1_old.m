function vdp() 
u  = 5; 
x0 = 1; 
y0 = 0; 
options=odeset();
[t,x] = ode45( @rhs, [0,25], [x0,y0],options); 
 
subplot(2,1,1); 
plot(t,x(:,1)); 
xlabel('t'); ylabel('x(t)'); 
title(sprintf('Van Der Pol, position vs, time, u=%3.2f',u)); 
 
%subplot(2,1,2); 
%plot(x(:,1),x(:,2)); 
%xlabel('x(t)'); ylabel('v(t)'); 
%hold on; 
%plot(x0,y0,'*r','MarkerSize',10); 
%title(sprintf('Phase portait showing limit cycle. x(0)=%3.2f, y(0)=%3.2f',x0,y0)); 
%axis equal 
%hold off; 
 
    function dxdt=rhs(t,x) 
        dxdt_1 = x(2); 
        dxdt_2 = u*(1-x(1)^2)*x(2)-x(1); 
 
        dxdt = [dxdt_1 ; dxdt_2]; 
    end 
end

% Defining the function:
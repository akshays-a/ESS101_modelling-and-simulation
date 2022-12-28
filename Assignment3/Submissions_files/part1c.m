x=sym('x',[3 ,1],'real');
syms u tk t  dt real
alpha=sym('alpha',[2,1],'real');
x=@(t)(x)
f=@(t)(x+alpha(1)*(t-tk)+alpha(2)*(t-tk)^2);

xdot=jacobian(f(tk),t)
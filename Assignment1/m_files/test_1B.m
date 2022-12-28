syms m1 m2 L theta phi thetaDot phiDot q2Dot g real
syms p q P1Dot P2Dot p11Dot p12Dot p13Dot p21Dot p22Dot p23Dot Z real

%Position of the helicopter is as follows
P1 = sym('p1%d',[1,3],'real')';

%First-order derivate of the position vector of helicopter
P1Dot = [p11Dot ; p12Dot ; p13Dot]; 

%Position of the hovering mass is as follows
P2= sym('p2%d',[1,3],'real')';

%First-order derivate of the position vector of hovering mass
P2Dot = [p21Dot ; p22Dot ; p23Dot];

%Generalized co-ordinates for the system are defined as following
q = [P1; P2];

%First-order derivate of the generailzed co-ordinate system
qDot = [P1Dot; P2Dot]; 

%Kinetic energy of the system is defined as T
T = simplify(0.5 * m1 * (P1Dot'* P1Dot) +  0.5 * m2 * (P2Dot'* P2Dot));

%Potential energy of the system is defined as V
V = m1*g*[0 0 1]*P1 + m2*g*[0 0 1]*P2;

%Constaints for the model are defined as following
e = P1 - P2;
C = 0.5 * (e' * e- L^2);
CDot = jacobian(C,q)*qDot;

%Lagrange equation is as follows
LF = simplify(T - V - Z*C);

%Derivative of Lagrange equation with qDot
dLdqDot = jacobian(LF, qDot);

%Derivative of Lagrange equation with q
dLdq = jacobian(LF, q);

%Time derivative of derivative of Lagrange equation with qDot
dt_dLdqDot = jacobian(dLdqDot,q)*qDot + jacobian(dLdqDot,qDot)*q2Dot;

%External forces acting on the helicopter are as follows
u = sym('u1%d',[1,3],'real')';

%Transpose of the derivative of position of helicopter with q
dP1dq_Transpose = jacobian(P1,q)';

%Amount of work produced on the system when moving in the generalized coordinates
Q = dP1dq_Transpose * u;

%When the equations are assembled for the Euler-Lagrange equation in the form: dt_dLdqDot - dLdq = Q
% and if v = qDot then the model can be written as M*vDot = b, where M and b are as follows
M = simplify(jacobian(dLdqDot,qDot));

b = simplify(Q - jacobian(dLdqDot,q)*qDot + dLdq);

%Solving for q2Dot and Z can be done by putting the model in the form 
% [M ,a(q) ; a(q)' , 0] [q2Dot ; Z] = c. Where a and c are as follows

a = jacobian(C,q);

c = simplify([Q-jacobian(dLdqDot,q)*qDot+jacobian(T,q)'-jacobian(V, q)'; - jacobian(CDot, q)*qDot]);

LHS = ([M, a; a', 0])\c;
LHS = simplify(LHS);

function y = pondcstr_OutputFcn(x, u)

% Initialize variables
h = x(1);
c = x(2);
theta = u(1);
q_in = u(2);
c_in = u(3);

co = 0.67;
Ao = 1.5;
g = 9.81;

y = theta*co*Ao*sqrt(2*g*h);
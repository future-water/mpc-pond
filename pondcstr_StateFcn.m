function dxdt = pondcstr_StateFcn(x, u)

% Initialize variables
h = x(1);
c = x(2);

theta = u(1);

q_in = u(2);
c_in = u(3);

dxdt = zeros(2,1);

% Parameters
co = 1;
Ao = 1;
g = 9.81;
A = 100;
k = 42.048/365/24;

% Model equations
%q_out = co*Ao*theta*sqrt(2*g*h);

dxdt(1) = 1/A*(q_in-co*Ao*theta*sqrt(2*g*h));
dxdt(2) = (-q_in/(A*h)-k)*c+q_in/(A*h)*c_in;




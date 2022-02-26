function x_dt = pondcstr_StateFcn(x, u)

% Initialize variables
h = x(1);
c = x(2);
theta = u(1);
q_in = u(2);
c_in = u(3);

x_dt = zeros(2,1);

% Parameters
co = 0.67;
Ao = 1.5;
g = 9.81;
A = 100;
k = 0.05;
dt = 1;

% Model equations
q_out = co*Ao*sqrt(2*g*h);

x_dt(1) = max( (h + dt/A*(q_in-theta*q_out)), 0);

if h > 0
     x_dt(2) = (c*A*h*exp(-k*dt)+c_in*q_in*dt)/(A*h+q_in*dt);
else
     x_dt(2) = 0;
end


function x_dt = pondcstr_StateFcn(x, u)
% Initialize variables
h = x(1);
c = x(2);

theta = u(1);

q_in = u(2);
c_in = u(3);

x_dt = zeros(2,1);

% Parameters
co = 0.65;
Ao = 1;
g = 32.2;
k = 0.9/24/60/60;
dt = 15*60;

% % node 4
h = max(0, h);
elevation = [0, 2, 4, 6, 8, 10];
area = [82971, 93258, 106100, 119152, 134285, 134285];
A = interp1(elevation, area, h,'spline');

% Model equations
q_out = theta*co*Ao*sqrt(2*g*h)*min(1,h);

x_dt(1) = max(0,h + dt/A*(q_in-q_out));

if h > 0 
    x_dt(2) = (c*A*h*exp(-k*dt)+c_in*q_in*dt)/(A*h+q_in*dt);
else
    x_dt(2) = 0;
end

end
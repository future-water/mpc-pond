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

% % node 2
% elevation = [0, 0.5, 1.5, 2.5, 2.7, 10];
% area = [0, 1360, 2046, 2960, 3166, 3166];
% A = interp1(elevation, area, h,'linear','extrap');

% % node 3
% elevation = [0, 2, 4, 6, 8, 10];
% area = [78572, 93640, 104101, 115196, 128463, 128463];
% A = interp1(elevation, area, h,'linear','extrap');

% % node 1
% elevation = [0, 1, 3, 5];
% area = [63669, 81169, 85949, 102846.16];
% A = interp1(elevation, area, h,'linear');

% node 4
elevation = [0, 2, 4, 6, 8, 10];
area = [82971, 93258, 106100, 119152, 134285, 134285];
A = interp1(elevation, area, h,'spline');

% Model equations
q_out = theta*co*Ao*sqrt(2*g*h)*min(1,h);

x_dt(1) = max(0,h + dt/A*(q_in-q_out));
%x_dt(2) = q_in/(A*h)*(c_in-c)-k*c;
%x_dt(2) = (c*A*h*exp(-k*dt)+c_in*q_in*dt)/(A*h+q_in*dt);

if h > 0 
    x_dt(2) = (c*A*h*exp(-k*dt)+c_in*q_in*dt)/(A*h+q_in*dt);
else
    x_dt(2) = 0;
end

end
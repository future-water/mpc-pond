function y = pondcstr_OutputFcn(x, u)

% Initialize variables
h = x(1);
c = x(2);

theta = u(1);

q_in = u(2);
c_in = u(3);

co = 0.65;
Ao = 1;
g = 32.2;

y = zeros(3,1);
y(1) = h;
y(2) = c;
y(3) = theta*co*Ao*sqrt(2*g*h)*min(1,h);
function y = pondcstr_OutputFcn(x, u)

% Initialize variables
h = x(1);
c = x(2);

theta = u(1);

q_in = u(2);
c_in = u(3);

y = theta*sqrt(2*9.81*h);
%% Initializing storm event runoff
load('Stormevent_Runoff.mat')
%%  Isolating the strom for analysis 
for i=1:800
    Temp1(i)=0;
end
Runoff=horzcat(transpose(A(1210:1737)),Temp1); % Extending the duration of dry period after storm in order to visualize the behavior of system after the event.
Flow=Runoff*0.0004719474*100;
%% Generating Pollutant Loads 
for i=1:length(Runoff)
    if Runoff(i)==0
        Nitrate(i)=0;
    else 
        Nitrate(i)=( Runoff(i)*0.0004719474*100 + 3);
    end
end
qin = Flow(1:801)';
cin = Nitrate(1:801)';


co = 1;
Ao = 1;
g = 9.81;
A = 100;
k = 42.048/365/24/60;
for i =1:length(qin)

    syms h(t) c(t)

ode1 = diff(h) == 1/A*(qin(i)-co*Ao*sqrt(2*g*h));
ode2 = diff(c) == (-qin(i)/(A*h)-k)*c+qin(i)/(A*h)*cin(i);
odes = [ode1; ode2];

S(i) = dsolve(odes);

end

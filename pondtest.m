clear all
clc
%% Generating Pollutant Loads 
Runoff = [0	0	0	0.148394421	0.221659535	0.341659535	1.278727229	2.766291368	2.631526788	3.339896238	4.787321158	5.610194486	6.446603266	7.230144498	8.000442886	7.397921799	6.446603266	5.832476992	4.303919595	4.282568694	3.686834211	3.046184493	2.86963369	2.565435274	2.581877922	1.952045239	1.853582852	1.798501871	1.599646834	1.437131747	1.388960076	1	0.5	0.3	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0.148394421	0.221659535	0.341659535	1.278727229	2.766291368	2.631526788	3.339896238	4.787321158	5.610194486	7.446603266	8.230144498	10.00044289	8.397921799	7.446603266	5.832476992	4.303919595	4.282568694	3.686834211	3.046184493	2.86963369	2.565435274	2.581877922	1.952045239	1.853582852	1.798501871	1.599646834	1.437131747	1.388960076	0.5	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
];

for i=1:length(Runoff)
    if Runoff(i)==0
        Nitrate(i)=0;
    else 
        Nitrate(i)=( Runoff(i)*0.0004719474*100 + 1);
%        Nitrate(i)=( Runoff(i)*0.5);
    end
end
qins = Runoff;
cins = Nitrate;

qins = [qins zeros(1,200)];
cins = [cins zeros(1,200)];

MD = [qins' cins'];

%% Nonlinear MPC Design 
% Create a nonlinear MPC object with 2 states, 1 outputs, 1 manipulated
% variables, and 1 measured disturbance.
nlmpcobj_Plan = nlmpc(2, 1, 'MV', 1, 'MD', [2,3]);

% Controller sample time |Ts| and prediction horizon.
Ts = 1;
nlmpcobj_Plan.Ts = Ts;
nlmpcobj_Plan.PredictionHorizon = 400;
nlmpcobj_Plan.ControlHorizon = [4*ones(1,100)];

% Specify the nonlinear model in the controller
nlmpcobj_Plan.Model.StateFcn = @(x,u) pondcstr_StateFcn(x, u);
nlmpcobj_Plan.Model.OutputFcn = @(x,u) pondcstr_OutputFcn(x,u);

% Initial condition
x0(1) = 0.001;
x0(2) = 0.001;
u0 = 1;

% To simulate no-control case, valve is always open
nlmpcobj_Plan.MV(1).Min = 1;
nlmpcobj_Plan.MV(1).Max = 1;
 
fprintf('\nOptimization started...\n');
[~,~,nc] = nlmpcmove(nlmpcobj_Plan,x0,u0,zeros(1,1),MD);
fprintf('Optimization finished...\n');

% control 
% Specify the bounds for orifice opening ratio 
nlmpcobj_Plan.MV(1).Min = 0;
nlmpcobj_Plan.MV(1).Max = 1;

% Specify the upper bound for the height of the pond
hlimit = 1;
nlmpcobj_Plan.State(1).Max = hlimit;
nlmpcobj_Plan.State(1).Min = 0;

% cost function
nlmpcobj_Plan.Optimization.CustomCostFcn = 'pondcstrCostFcn'; 
nlmpcobj_Plan.Optimization.ReplaceStandardCost = true;

%To configure the manipulated variables to vary linearly with time
nlmpcobj_Plan.Optimization.MVInterpolationOrder = 1;

% Find the optimal trajectory for the manipulated variable such that the total pollution load is minimized 
fprintf('\nOptimization started...\n');
[~,~,Info] = nlmpcmove(nlmpcobj_Plan,x0,u0,zeros(1,1),MD);
fprintf('   First order optimality (Info.ExitFlag = %i).\n',Info.ExitFlag);
fprintf('Optimization finished...\n');

% plotting
Nstep = size(Info.Xopt,1) - 1;
figure 
subplot(2,1,1)
plot((1:length(qins))*Ts,qins)
title('q_in')
subplot(2,1,2)
plot((1:Nstep)*Ts, cins)
title('c_in')

figure
subplot(3,2,1)
plot((0:Nstep)*Ts,nc.Xopt(:,1),(0:Nstep)*Ts, Info.Xopt(:,1))
legend({'No control','Control'})
title('height (y1)')
subplot(3,2,2)
plot((0:Nstep)*Ts, nc.Xopt(:,2),(0:Nstep)*Ts, Info.Xopt(:,2))
legend({'No control','Control'})
title('concentration (y2)')

subplot(3,2,3)
plot((0:Nstep)*Ts,nc.Yopt(:,1),(0:Nstep)*Ts,Info.Yopt(:,1))
legend({'No control','Control'})
title('outflow')
subplot(3,2,4)
plot((0:Nstep)*Ts, nc.MVopt(:,1),(0:Nstep)*Ts, Info.MVopt(:,1))
title('valve opening')

subplot(3,2,5)
plot((0:Nstep)*Ts, cumsum(nc.Xopt(:,2).*nc.Yopt(:,1)),(0:Nstep)*Ts, cumsum(Info.Xopt(:,2).*Info.Yopt(:,1)))
legend({'No control','Control'})
title('total pollutant loads')

subplot(3,2,6)
plot((0:Nstep)*Ts,cumsum((Info.Yopt(:,1)-mean(Info.Yopt(:,1))).^2),(0:Nstep)*Ts,cumsum((Info.Xopt(:,2).*Info.Yopt(:,1)).^2),(0:Nstep)*Ts, cumsum((Info.Yopt(:,1)-mean(Info.Yopt(:,1))).^2+(Info.Xopt(:,2).*Info.Yopt(:,1)).^2))
legend({'minimizing peak flow','minimizing pollutant loads','Total'})
title('Cumulative Costs')


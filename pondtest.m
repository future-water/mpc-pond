clear
clc
%% Nonlinear MPC Design 
% Create a nonlinear MPC object with 2 states, 1 outputs, 1 manipulated
% variables, and 1 measured disturbance.
nlmpcobj_Plan = nlmpc(2, 1, 'MV', 1, 'MD', [2,3]);

% Given the expected batch duration |Tf|, choose the controller sample time
% |Ts| and prediction horizon.
Ts = 1;
nlmpcobj_Plan.Ts = Ts;
nlmpcobj_Plan.PredictionHorizon = 200;
nlmpcobj_Plan.ControlHorizon = 30;

%nlmpcobj_Plan.Model.StateFcn = @(x,u) pondcstr_StateFcnDT(x,u,Ts);
nlmpcobj_Plan.Model.StateFcn = @(x,u) pondcstr_StateFcn(x, u);
nlmpcobj_Plan.Model.OutputFcn = @(x,u) pondcstr_OutputFcn(x,u);
%nlmpcobj_Plan.Model.IsContinuousTime = false;

% Specify the bounds for orifice opening ratio
nlmpcobj_Plan.MV(1).Min = 0;
nlmpcobj_Plan.MV(1).Max = 1;

% Specify the bounds for the height of the pond
hlimit = 1.5;
nlmpcobj_Plan.State(1).Max = hlimit;
%nlmpcobj_Plan.State(2).Min = 0;
nlmpcobj_Plan.OV.Max = 2;
%nlmpcobj_Plan.OV(1).Min = 0;

% cost function
%nlmpcobj_Plan.Optimization.CustomCostFcn = 'pondcstrCostFcn'; 
%nlobj_Plan.Optimization.CustomCostFcn = 'pondcstrCostFcn';
%nlmpcobj_Plan.Optimization.ReplaceStandardCost = true;

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

for i = 1:100
     qins(i) = qin(8*i);
     cins(i) = cin(8*i);
end

qins = [qins zeros(1,100)];
cins = [cins zeros(1,100)];

MD = [qins' cins'];
load('nocontrol.mat')

x0(1) = 0.1;
x0(2) = 0.001;

u0 = 1;

fprintf('\nOptimization started...\n');
[~,~,Info] = nlmpcmove(nlmpcobj_Plan,x0,u0,zeros(1,1),MD);
fprintf('   Expected production of C (y1) is %g moles.\n',Info.Yopt(end,1));
fprintf('   First order optimality is satisfied (Info.ExitFlag = %i).\n',...
    Info.ExitFlag);
fprintf('Optimization finished...\n');
nlmpcobj_Plan.MV(1).Min = 1;
nlmpcobj_Plan.MV(1).Max = 1;
[~,~,nc] = nlmpcmove(nlmpcobj_Plan,x0,u0,zeros(1,1),MD);

Nstep = size(Info.Xopt,1) - 1;
figure 
subplot(2,1,1)
plot((1:length(qins))*Ts,qins)
title('q_in')
subplot(2,1,2)
plot((1:Nstep)*Ts, cins)
title('c_in')

figure
subplot(2,1,1)
%plot((0:Nstep)*Ts, Info.Yopt(:,1))
plot((0:length(qins))*Ts,h_nc,(0:Nstep)*Ts, Info.Xopt(:,1))
legend({'No control','Control'})
title('height (y1)')
subplot(2,1,2)
tTs = (0:Nstep)*Ts;
plot((0:length(qins))*Ts,c_nc,(0:Nstep)*Ts, Info.Xopt(:,2))
legend({'No control','Control'})
title('concentration (y2)')

figure
subplot(2,1,1)
plot((0:length(qins))*Ts,qout_nc,(0:Nstep)*Ts,Info.Yopt(:,1))
legend({'No control','Control'})
title('outflow')
subplot(2,1,2)
plot((0:Nstep)*Ts, Info.MVopt(:,1), (0:Nstep)*Ts, Info.Yopt(:,1)./sqrt(2*9.81*Info.Xopt(:,1)))
title('valve opening')

figure 
plot((0:length(qins))*Ts,cumsum(qout_nc.*c_nc),(0:Nstep)*Ts, cumsum(Info.Xopt(:,2).*Info.Yopt(:,1)))
legend({'No control','Control'})
title('total loads')

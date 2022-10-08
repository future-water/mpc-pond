clear all
close all
clc
%% Data preparation
% Importing input data
data1 = readtable('gamma_4_15min_base.csv'); % true
data2 = readtable('gamma_4_15min_200b_50rate_150emc.csv'); % forecasts

% Extending the duration to see the behavior of system after the storm event
qin_t = [data1.qin; zeros(500,1)]; 
cin_t = [data1.cin; zeros(500,1)];

qin_f = [data2.qin; zeros(500,1)];
%cin_f = [data2.cin; zeros(500,1)]; 
cin_f25 = [0.25*sum(data1.qin.*data1.cin)/sum(data1.qin)*ones(1077,1)]; % Imperfect water quality prediction as EMC
cin_f50 = [0.5*sum(data1.qin.*data1.cin)/sum(data1.qin)*ones(1077,1)]; 
cin_f75 = [0.75*sum(data1.qin.*data1.cin)/sum(data1.qin)*ones(1077,1)]; 
cin_f100 = [1*sum(data1.qin.*data1.cin)/sum(data1.qin)*ones(1077,1)]; 
cin_f125 = [1.25*sum(data1.qin.*data1.cin)/sum(data1.qin)*ones(1077,1)]; 
cin_f150 = [1.5*sum(data1.qin.*data1.cin)/sum(data1.qin)*ones(1077,1)]; 
cin_f175 = [1.75*sum(data1.qin.*data1.cin)/sum(data1.qin)*ones(1077,1)]; 
cin_f200 = [2*sum(data1.qin.*data1.cin)/sum(data1.qin)*ones(1077,1)]; 

MD_t = [qin_t, cin_t]; % Measured disturbances with perfect knowledge
cin_f = [cin_f25, cin_f50, cin_f75, cin_f100, cin_f125, cin_f150, cin_f175, cin_f200];

%% Nonlinear MPC Design 
% Create a nonlinear MPC object with 2 states, 3 outputs, 1 manipulated variables, and 2 measured disturbance.
nlmpcobj_Plan = nlmpc(2, 3, 'MV', 1, 'MD', [2,3]);

% Controller sample time |Ts| and prediction horizon.
Ts = 1;
nlmpcobj_Plan.Ts = Ts;
nlmpcobj_Plan.PredictionHorizon = length(MD_t);
nlmpcobj_Plan.ControlHorizon = 2;

% Specify the nonlinear model in the controller
nlmpcobj_Plan.Model.StateFcn = @(x,u) pondcstr_StateFcn(x, u);
nlmpcobj_Plan.Model.IsContinuousTime = false;
nlmpcobj_Plan.Model.OutputFcn = @(x,u) pondcstr_OutputFcn(x,u);

% Initial condition
x0(1) = 0.01;
x0(2) = 0;
u0 = 1;

%% Passive system simulation; valve is always opened
nlmpcobj_Plan.MV(1).Min = 1;
nlmpcobj_Plan.MV(1).Max = 1;

fprintf('\nPassive system\n');
[~,~,nc] = nlmpcmove(nlmpcobj_Plan,x0,u0,[],MD_t);
fprintf('Passive system End\n');

%% MPC simulation
horizon = 96; % 24hr
nlmpcobj_Plan.PredictionHorizon = horizon;
nlmpcobj_Plan.ControlHorizon = 2;

% Specify the bounds for orifice opening ratio 
nlmpcobj_Plan.MV(1).Min = 0;
nlmpcobj_Plan.MV(1).Max = 1;

% Specify the upper bound for the height of the pond
hlimit = 10;
nlmpcobj_Plan.State(1).Max = hlimit;

nlmpcobj_Plan.State(1).ScaleFactor = hlimit;
nlmpcobj_Plan.State(2).ScaleFactor = 25;

% Cost function
nlmpcobj_Plan.Optimization.CustomCostFcn = 'pondcstrCostFcn'; 
nlmpcobj_Plan.Optimization.ReplaceStandardCost = true;
yref = [0 0 0];

% Extended Kalman Filter
DStateFcn = @(xk,uk) pondcstr_StateFcn(xk, uk);
DMeasFcn = @(xk) xk(2);
EKF = extendedKalmanFilter(DStateFcn,DMeasFcn,x0);
EKF.ProcessNoise = diag([0;1]);
EKF.MeasurementNoise = 0.1;

%% MPC-false
% Find the optimal trajectory for valve opening
fprintf('\nMPC Optimization started...\n');
tic

falhh = [];
falcc = [];
falyy = [];
falmvmv = [];

for con = 1:8
    con
    % initial values
    falh(1) = x0(1);
    falc(1) = x0(2);
    faly(1) = 0;
    falmv(1) = u0;

waitbar_h = waitbar(0,'Process . . . ');
for k = 1:(length(MD_t)-horizon)
    % MPC-flase
    waitbar(k/(length(MD_t)-horizon),waitbar_h)
    [~,~,false] = nlmpcmove(nlmpcobj_Plan,[falh(k), falc(k)],falmv(k),yref,[qin_f(k:(horizon+k-1),:), cin_f(k:(horizon+k-1),con)]);
    falh(k+1,1) = false.Xopt(2,1);
    falc(k+1,1) = false.Xopt(2,2);
    faly(k+1,1) = false.Yopt(2,3);
    falmv(k+1,1) = false.MVopt(2,:);
end
close(waitbar_h); clear waitbar_h;
fprintf('MPC Optimization finished...\n');
timeElapsed = toc

falhh = [falhh, falh];
falcc = [falcc, falc];
falyy = [falyy, faly];
falmvmv = [falmvmv, falmv];
end

%% MPC-true
% Find the optimal trajectory for valve opening
fprintf('\nMPC Optimization started...\n');
tic
% initial values
trueh(1) = x0(1);
truec(1) = x0(2);
truey(1) = 0;
truemv(1) = u0;
waitbar_h = waitbar(0,'Process . . . ');
for k = 1:(length(MD_t)-horizon)
    waitbar(k/(length(MD_t)-horizon),waitbar_h)
    % MPC-true
    [~,~,true] = nlmpcmove(nlmpcobj_Plan,[trueh(k), truec(k)],truemv(k),yref,MD_t(k:(horizon+k-1),:));
    trueh(k+1,1) = true.Xopt(2,1);
    truec(k+1,1) = true.Xopt(2,2);
    truey(k+1,1) = true.Yopt(2,3);
    truemv(k+1,1) = true.MVopt(2,:);
end
close(waitbar_h); clear waitbar_h;
fprintf('MPC Optimization finished...\n');
timeElapsed = toc

%% MPC-EKF
% Find the optimal trajectory for valve opening with MPC-EKF
fprintf('\nMPC-EKF Optimization started...\n');
tic
hh = [];
cc = [];
yy = [];
mvmv = [];
faltruecc = [];

for con = 1:8
[~,~,mpcekf] = nlmpcmove(nlmpcobj_Plan,x0,u0,yref,[qin_f(1:horizon,:), cin_f(1:horizon,con)]);
mpc_ekfh = mpcekf.Xopt(1:2,1);
mpc_ekfc = mpcekf.Xopt(1:2,2);
mpc_ekfy = mpcekf.Yopt(1:2,3);
mpc_ekfmv = mpcekf.MVopt(1:2,:);
mpc_x = x0;
mpc_x(2,:) = pondcstr_StateFcn(mpc_x(1,:), [mpcekf.MVopt(1,:), MD_t(1,:)]);

EKF.ProcessNoise = diag([0;1]);
EKF.MeasurementNoise = 0.1;

waitbar_h = waitbar(0,'Process . . . ');
for k = 1:(length(MD_t)-horizon-1)
    waitbar(k/(length(MD_t)-horizon-1),waitbar_h)
    yk = pondcstr_StateFcn(mpc_x(k+1,:), [mpcekf.MVopt(2,:), MD_t(k+1,:)]);
    yk = max(0, yk);
    xxk = correct(EKF, yk(2));
    xxk = [mpcekf.Xopt(2,1), xxk(2)];
    [~,~,mpcekf] = nlmpcmove(nlmpcobj_Plan,xxk, mpcekf.MVopt(2,:),yref,[qin_f(k+1:(horizon+k),:), cin_f(k+1:(horizon+k),con)]);
    predict(EKF,[mpcekf.MVopt(2,:), qin_f(k+2,:), cin_f(k+2,con)]);

    mpc_ekfh = [mpc_ekfh; mpcekf.Xopt(2,1)];
    mpc_ekfc = [mpc_ekfc; mpcekf.Xopt(2,2)];
    mpc_ekfy = [mpc_ekfy; mpcekf.Yopt(2,3)];
    mpc_ekfmv = [mpc_ekfmv; mpcekf.MVopt(2,:)];
    mpc_x = [mpc_x; xxk];
end
close(waitbar_h); clear waitbar_h;
hh = [hh, mpc_ekfh];
cc = [cc, mpc_ekfc];
yy = [yy, mpc_ekfy];
mvmv = [mvmv, mpc_ekfmv];
faltruecc = [faltruecc, mpc_x(:,2)];

end

fprintf('MPC-EKF Optimization finished...\n');
timeElapsed = toc

%% Unit conversion from US to SI
qin_t = qin_t*0.028316846592;

nc.Xopt(:,1) = nc.Xopt(:,1)/3.281;
nc.Yopt(:,3) = nc.Yopt(:,3)*0.028316846592;

hlimit = hlimit/3.281;

falh = falh/3.281;
faly = faly*0.028316846592;

falhh = falhh/3.281;
falyy = falyy*0.028316846592;

trueh = trueh/3.281;
truey = truey*0.028316846592;

mpc_ekfh = mpc_ekfh/3.281;
mpc_ekfy = mpc_ekfy*0.028316846592;

hh = hh/3.281;
yy = yy*0.028316846592;

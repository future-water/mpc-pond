clear all
clc
tic
%% Generating Pollutant Loads

data1 = readtable('gamma_4_15min_x1.csv'); % true
data2 = readtable('gamma_4_15min_x0.1.csv'); % forecasts 

qin_t = [data1.qin; zeros(500,1)];
cin_t = [data1.cin; zeros(500,1)];

qin_f = [data2.qin];
cin_f = [data2.cin];

MD = [qin_f, cin_f];
MD = [MD; zeros(500,2)]*1;
%% Nonlinear MPC Design 
% Create a nonlinear MPC object with 2 states, 1 outputs, 1 manipulated
% variables, and 1 measured disturbance.
nlmpcobj_Plan = nlmpc(2, 3, 'MV', 1, 'MD', [2,3]);

% Controller sample time |Ts| and prediction horizon.
Ts = 1;
nlmpcobj_Plan.Ts = Ts;
nlmpcobj_Plan.PredictionHorizon = length(MD);
nlmpcobj_Plan.ControlHorizon = 2;

% Specify the nonlinear model in the controller
nlmpcobj_Plan.Model.StateFcn = @(x,u) pondcstr_StateFcn(x, u);
nlmpcobj_Plan.Model.IsContinuousTime = false;
nlmpcobj_Plan.Model.OutputFcn = @(x,u) pondcstr_OutputFcn(x,u);

% Initial condition
x0(1) = 0.01;
x0(2) = 0;
u0 = 1;

% To simulate no-control case, valve is always open
nlmpcobj_Plan.MV(1).Min = 1;
nlmpcobj_Plan.MV(1).Max = 1;

fprintf('\nOptimization started...\n');
[~,~,nc] = nlmpcmove(nlmpcobj_Plan,x0,u0,[],MD);
fprintf('Optimization finished...\n');

%control
horizon = 96;
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

% cost function
nlmpcobj_Plan.Optimization.CustomCostFcn = 'pondcstrCostFcn'; 
nlmpcobj_Plan.Optimization.ReplaceStandardCost = true;

%nlmpcobj_Plan.Optimization.CustomIneqConFcn = @(X,U,e,data) [X((end),1) - 0.1];

%To configure the manipulated variables to vary linearly with time
%nlmpcobj_Plan.Optimization.MVInterpolationOrder = 1;
nlmpcobj_Plan.Optimization.UseSuboptimalSolution = true;
yref = [0, 0, 0];
u0 = 1;

% Kalman filter
DStateFcn = @(xk,uk) pondcstr_StateFcn(xk, uk);
%DMeasFcn = @(xk) xk(1:3);
DMeasFcn = @(xk) xk(2);

EKF = extendedKalmanFilter(DStateFcn,DMeasFcn,x0);

%EKF.ProcessNoise = diag([1e-6;0.1;0]);
%EKF.MeasurementNoise = diag([0.01;0.01;0]);
% EKF.ProcessNoise = diag([0;1]);
% EKF.MeasurementNoise = diag([0.01;0.01]);
EKF.ProcessNoise = 1;
EKF.MeasurementNoise = 0.01;

% Find the optimal trajectory for the manipulated variable such that the total pollution load is minimized 
fprintf('\nOptimization started...\n');
[~,~,Info] = nlmpcmove(nlmpcobj_Plan,x0,u0,yref,MD(1:horizon,:));
[~,~,true] = nlmpcmove(nlmpcobj_Plan,x0,u0,yref,[qin_t(1:horizon), cin_t(1:horizon)]);
fprintf('   First order optimality (Info.ExitFlag = %i).\n',Info.ExitFlag);
fprintf('Optimization finished...\n');
timeElapsed = toc

%predicted 
orinh = Info.Xopt(1:2,1);
orinc = Info.Xopt(1:2,2);
oriny = Info.Yopt(1:2,3);
orinmv = Info.MVopt(1:2,:);


trueh = true.Xopt(1:2,1);
truec = true.Xopt(1:2,2);
truey = true.Yopt(1:2,3);
truemv = true.MVopt(1:2,:);
randomness = [];

% %measurements
% yk = pondcstr_StateFcn(x0, [Info.MVopt(1,:), qin_t(1), cin_t(1)])';
% xk = correct(EKF, yk);

for k = 1:(length(MD)-horizon-1)
    % one control
    fprintf('\nOptimization started...\n');
    random = [0, 0.01].*randn(1,2);
    next = max(0, true.Xopt(2,:) + random);
    [~,~,Info] = nlmpcmove(nlmpcobj_Plan,[max(0,Info.Xopt(2,1)+0*randn(1)), max(0,Info.Xopt(2,2)+0*randn(1))], Info.MVopt(2,:),yref,MD(k+1:(horizon+k),:));
    [~,~,true] = nlmpcmove(nlmpcobj_Plan,next, true.MVopt(2,:),yref,[qin_t(k+1:(horizon+k),:), cin_t(k+1:(horizon+k),:)]);
    fprintf('   First order optimality (Info.ExitFlag = %i).\n',Info.ExitFlag);
    fprintf('Optimization finished...\n');
    k
    %openx = [openx, Info.Xopt];
    orinh = [orinh; Info.Xopt(2,1)];
    orinc = [orinc; Info.Xopt(2,2)];
    oriny = [oriny; Info.Yopt(2,3)];
    orinmv = [orinmv; Info.MVopt(2,:)];

    trueh = [trueh; true.Xopt(2,1)];
    truec = [truec; true.Xopt(2,2)];
    truey = [truey; true.Yopt(2,3)];
    truemv = [truemv; true.MVopt(2,:)];
    randomness = [randomness; random];

end


hh = [];
cc = [];
yy = [];
mvmv = [];
randrand = [];

for t = 1:1

% Find the optimal trajectory for the manipulated variable such that the total pollution load is minimized 
fprintf('\nOptimization started...\n');
[~,~,Info] = nlmpcmove(nlmpcobj_Plan,x0,u0,yref,MD(1:horizon,:));
fprintf('   First order optimality (Info.ExitFlag = %i).\n',Info.ExitFlag);
fprintf('Optimization finished...\n');
timeElapsed = toc

openx = Info.Xopt;
openy = Info.Yopt;
openmv = Info.MVopt;

inith = Info.Xopt(1:2,1);
initc = Info.Xopt(1:2,2);
inity = Info.Yopt(1:2,3);
initmv = Info.MVopt(1:2,:);
xx = [x0];

for k = 1:(length(MD)-horizon-1)
    % one control
    fprintf('\nOptimization started...\n');
    yk = pondcstr_StateFcn(xx(k,:), [Info.MVopt(2,:), qin_t(k), cin_t(k)]);
    yk = max(0, yk+randomness(k,:));
    xk = correct(EKF, yk(2));
    xk = [Info.Xopt(2,1), xk(2)];
    %error = xk - yk;
    [~,~,Info] = nlmpcmove(nlmpcobj_Plan,xk, Info.MVopt(2,:),yref,MD(k+1:(horizon+k),:));
    %[~,~,Info] = nlmpcmove(nlmpcobj_Plan,[max(0,Info.Xopt(2,1)+0*randn(1)), max(0,Info.Xopt(2,2)+0*randn(1))], Info.MVopt(2,:),yref(k+1:(horizon+k),:),MD(k+1:(horizon+k),:));
    fprintf('   First order optimality (Info.ExitFlag = %i).\n',Info.ExitFlag);
    fprintf('Optimization finished...\n');
    predict(EKF,[Info.MVopt(2,:),MD(k+2,:)]);
    k
    %openx = [openx, Info.Xopt];
    inith = [inith; Info.Xopt(2,1)];
    initc = [initc; Info.Xopt(2,2)];
    inity = [inity; Info.Yopt(2,3)];
    initmv = [initmv; Info.MVopt(2,:)];
    xx = [xx;xk];
    %randomness = [randomness; random];
    %er = [er; error];
end

hh = [hh, inith];
cc = [cc, initc];
yy = [yy, inity];
mvmv = [mvmv, initmv];
end


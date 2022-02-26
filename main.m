clear all
clc

% hydrograph
qins = [0 1 2 3 4 5 6 7 8 9 10 9 8 7 6 5 4 3 2 1];
qins = [qins zeros(1,80)];
%pollutograph
cins = 2*qins;
MD = [qins' cins'];

% Trajectory Planning
% Create a nonlinear MPC object with 2 states, 1 outputs, 1 manipulated
% variables, and 1 measured disturbance.
nlmpcobj_Plan = nlmpc(2, 1, 'MV', 1, 'MD', [2,3]);

% Controller sample time |Ts| and prediction horizon.
Ts = 1;
nlmpcobj_Plan.Ts = Ts;
nlmpcobj_Plan.PredictionHorizon = 100;
nlmpcobj_Plan.ControlHorizon = 100;

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

% mpc
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
nlmpcobj_Plan.State(2).Min = 0;

nlmpcobj_Plan.State(1).ScaleFactor = 1;
nlmpcobj_Plan.State(2).ScaleFactor = 15;


% cost function
nlmpcobj_Plan.Optimization.CustomCostFcn = 'pondcstrCostFcn'; 
nlmpcobj_Plan.Optimization.ReplaceStandardCost = true;

% pond will be fully drained at the end
nlmpcobj_Plan.Optimization.CustomIneqConFcn = @(X,U,e,data) [X(end,1)' - 0.0001];

%To configure the manipulated variables to vary linearly with time
nlmpcobj_Plan.Optimization.MVInterpolationOrder = 1;

% Find the optimal trajectory for the manipulated variable such that the total pollution load is minimized 
fprintf('\nOptimization started...\n');
[~,~,Info] = nlmpcmove(nlmpcobj_Plan,x0,u0,zeros(1,1),MD);
fprintf('   First order optimality (Info.ExitFlag = %i).\n',Info.ExitFlag);
fprintf('Optimization finished...\n');


%% reculsive 

EKF = extendedKalmanFilter(@pondcstr_StateFcn,@pondcstr_MeasurementFcn, x0);
EKF.MeasurementNoise = 0.01;

Tsteps = 100;        
xHistory = x0;
uHistory = [];
lastMV = zeros(1,1);

p =100;

options = nlmpcmoveopt;
for k = 1:Tsteps
    % Obtain plant output measurements with sensor noise.
    %yk = xHistory(k,1:2)' + randn*0.01;
    yk = xHistory(k,1:2)' ;
    % Correct state estimation based on the measurements.
    %xk = correct(EKF, yk);
    xk = yk;
    % Compute the control moves with reference previewing.
    [uk,options] = nlmpcmove(nlmpcobj_Plan,xk,lastMV,[],MD(k:min(k+9,Tsteps),:),options);
    % Predict the state for the next step.
    %predict(EKF,uk);
    % Store the control move and update the last MV for the next step.
    uHistory(k,:) = uk';
    lastMV = uk;
    
    x = pondcstr_StateFcn(xk, [uk MD(k, :)]);

    % Save plant states
    xHistory = [xHistory; x'];
end




clear all
close all
clc
%% Data preparation
% Importing input data
%data1 = readtable('gamma_4_15min_base.csv'); % true
data1 = readtable('pond4_background.csv'); % true
data2 = readtable('pond4_background.csv'); % forecasts

% Extending the duration to see the behavior of system after the storm event
qin_t = [data1.qin; zeros(789,1)]; 
cin_t = [data1.cin; zeros(789,1)];

qin_f = [data2.qin; zeros(789,1)];
cin_f = [sum(data1.qin.*data1.cin)/sum(data1.qin)*ones(1077,1)]; % Imperfect water quality prediction as EMC

MD_t = [qin_t, cin_t]; % Measured disturbances with perfect knowledge
MD = [qin_f, cin_f]; % Measured disturbances with imperfect knowledge

%% Nonlinear MPC Design 
% Create a nonlinear MPC object with 2 states, 3 outputs, 1 manipulated variables, and 2 measured disturbance.
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

%% MPC-true and MPC-false
% Find the optimal trajectory for valve opening
fprintf('\nMPC Optimization started...\n');
tic
% initial values
falh(1) = x0(1);
falc(1) = x0(2);
faly(1) = 0;
falmv(1) = u0;

trueh(1) = x0(1);
truec(1) = x0(2);
truey(1) = 0;
truemv(1) = u0;

waitbar_h = waitbar(0,'Process . . . ');
for k = 1:(length(MD)-horizon)
    % MPC-flase
    waitbar(k/(length(MD)-horizon),waitbar_h)
    [~,~,false] = nlmpcmove(nlmpcobj_Plan,[falh(k), falc(k)],falmv(k),yref,MD(k:(horizon+k-1),:));
    falh(k+1,1) = false.Xopt(2,1);
    falc(k+1,1) = false.Xopt(2,2);
    faly(k+1,1) = false.Yopt(2,3);
    falmv(k+1,1) = false.MVopt(2,:);

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
[~,~,mpcekf] = nlmpcmove(nlmpcobj_Plan,x0,u0,yref,MD(1:horizon,:));
mpc_ekfh = mpcekf.Xopt(1:2,1);
mpc_ekfc = mpcekf.Xopt(1:2,2);
mpc_ekfy = mpcekf.Yopt(1:2,3);
mpc_ekfmv = mpcekf.MVopt(1:2,:);
mpc_x = x0;
mpc_x(2,:) = pondcstr_StateFcn(mpc_x(1,:), [mpcekf.MVopt(1,:), MD_t(1,:)]);

waitbar_h = waitbar(0,'Process . . . ');
for k = 1:(length(MD)-horizon-1)
    waitbar(k/(length(MD)-horizon-1),waitbar_h)
    yk = pondcstr_StateFcn(mpc_x(k+1,:), [mpcekf.MVopt(2,:), MD_t(k+1,:)]);
    xxk = correct(EKF, yk(2));
    xxk = [mpcekf.Xopt(2,1), xxk(2)];
    [~,~,mpcekf] = nlmpcmove(nlmpcobj_Plan,xxk, mpcekf.MVopt(2,:),yref,MD(k+1:(horizon+k),:));
    predict(EKF,[mpcekf.MVopt(2,:),MD(k+2,:)]);

    mpc_ekfh = [mpc_ekfh; mpcekf.Xopt(2,1)];
    mpc_ekfc = [mpc_ekfc; mpcekf.Xopt(2,2)];
    mpc_ekfy = [mpc_ekfy; mpcekf.Yopt(2,3)];
    mpc_ekfmv = [mpc_ekfmv; mpcekf.MVopt(2,:)];
    mpc_x = [mpc_x; xxk];
end
close(waitbar_h); clear waitbar_h;
fprintf('MPC-EKF Optimization finished...\n');
timeElapsed = toc

%% Multiple simulations with noise
hh = [];
cc = [];
yy = [];
mvmv = [];
randrand = [];

for t = 1:1 % the number of iterations for multiple sumulation with noise
[~,~,mpcekf] = nlmpcmove(nlmpcobj_Plan,x0,u0,yref,MD(1:horizon,:));
mpc_ekfh_n = mpcekf.Xopt(1:2,1);
mpc_ekfc_n = mpcekf.Xopt(1:2,2);
mpc_ekfy_n = mpcekf.Yopt(1:2,3);
mpc_ekfmv_n = mpcekf.MVopt(1:2,:);
mpc_x_n = x0;
mpc_x_n(2,:) = pondcstr_StateFcn(mpc_x_n(1,:), [mpcekf.MVopt(1,:), MD_t(1,:)]);
randomness = [];
measure_c = mpcekf.Xopt(1:2,2);

EKF.ProcessNoise = diag([0;1]);
EKF.MeasurementNoise = 0.1;

for k = 1:(length(MD)-horizon-1)
    yk_n = pondcstr_StateFcn(mpc_x_n(k+1,:), [mpcekf.MVopt(2,:), MD_t(k+1,:)]);
    random = [0, 0.02].*randn(1,2).*yk_n';
    yk_n = max(0, yk_n'+random);
    xxk = correct(EKF, yk_n(2));
    xxk = [mpcekf.Xopt(2,1), xxk(2)];
    [~,~,mpcekf] = nlmpcmove(nlmpcobj_Plan,xxk, mpcekf.MVopt(2,:),yref,MD(k+1:(horizon+k),:));
    predict(EKF,[mpcekf.MVopt(2,:),MD(k+2,:)]);

    mpc_ekfh_n = [mpc_ekfh_n; mpcekf.Xopt(2,1)];
    mpc_ekfc_n = [mpc_ekfc_n; mpcekf.Xopt(2,2)];
    measure_c = [measure_c; yk_n(2)];
    mpc_ekfy_n = [mpc_ekfy_n; mpcekf.Yopt(2,3)];
    mpc_ekfmv_n = [mpc_ekfmv_n; mpcekf.MVopt(2,:)];
    mpc_x_n = [mpc_x_n; xxk];
    randomness = [randomness; random];
end

hh = [hh, mpc_ekfh_n];
cc = [cc, mpc_ekfc_n];
yy = [yy, mpc_ekfy_n];
mvmv = [mvmv, mpc_ekfmv_n];
randrand = [randrand, randomness];
end

%% RBC-outflow
% Parameters
co = 0.65;
Ao = 1;
g = 32.2;
k = 0.8/24/60/60;
dt = 15*60;

h0 = 0.01; % S(0)
C0 = 0;
h_q(1) = h0;
C_q(1) = C0;
q_in = qin_t;
c_in = cin_t;

T = length(trueh);
q_desired = max(truey);
%q_desired = 2;

hdes = 1/(2*g)*(q_desired/(co*Ao))^2;
hlimit = max(trueh);
indicator = 1;
h_reten = 0.05;

% Operating policy
for t = 1:T
    % pond configuration
    elevation = [0, 2, 4, 6, 8, 10];
    area = [82971, 93258, 106100, 119152, 134285, 134285];
    A = interp1(elevation, area, h_q(t),'spline');

    C_q(t+1) = (C_q(t)*A*h_q(t)*exp(-k*dt)+c_in(t)*q_in(t)*dt)/(A*h_q(t)+q_in(t)*dt);

    if (h_q(t) < hlimit) && (indicator == 1)
        q_out_q(t) = 0;
        h_q(t+1) = h_q(t) + dt/A*(q_in(t) - q_out_q(t));
        theta_q(t) = 0;
        indicator = 1;

    elseif (h_q(t) >= hlimit)
        q_out_q(t) = q_desired;
        h_q(t+1) = max(0, h_q(t) + dt/A*(q_in(t) - q_out_q(t)));
        theta_q(t) = 1;
        indicator = 0;

    elseif (h_q(t) < hlimit) && (h_q(t) >= hdes) && (indicator == 0)
        q_out_q(t) = q_desired;
        h_q(t+1) = max(0, h_q(t) + dt/A*(q_in(t) - q_out_q(t)));
        theta_q(t) = q_desired/(co*Ao*sqrt(2*g*h_q(t)));
        indicator = 0;

    elseif (h_q(t) < hdes) && (indicator == 0) && (h_q(t) > h_reten)
        q_out_q(t) = co*Ao*sqrt(2*g*h_q(t))*min(1,h_q(t));
        h_q(t+1) = max(0, h_q(t) + dt/A*(q_in(t) - q_out_q(t)));
        theta_q(t) = 1;
        indicator = 0;

    elseif (h_q(t) <= h_reten) && (indicator == 0)
        q_out_q(t) = 0;
        h_q(t+1) = h_q(t) + dt/A*(q_in(t) - q_out_q(t));
        theta_q(t) = 0;
        indicator = 1;
    end
end
spill_q = 134285/10.764*max(0, h_q-10)/3.281/15/60.* C_q; 

%% RBC-Concentration

climit = 5;
h_c(1) = h0;
C_c(1) = C0;
% Operating policy
for t = 1:T
    % pond configuration
    elevation = [0, 2, 4, 6, 8, 10];
    area = [82971, 93258, 106100, 119152, 134285, 134285];
    A = interp1(elevation, area, h_c(t),'spline');

    C_c(t+1) = (C_c(t)*A*h_c(t)*exp(-k*dt)+c_in(t)*q_in(t)*dt)/(A*h_c(t)+q_in(t)*dt);

    if (C_c(t) > climit) && (h_c(t) < hlimit)
        q_out_c(t) = 0;
        h_c(t+1) = h_c(t) + dt/A*(q_in(t) - q_out_c(t));
        theta_c(t) = 0;

    elseif (h_c(t) >= hlimit)
        q_out_c(t) = co*Ao*sqrt(2*g*h_c(t))*min(1,h_c(t));
        h_c(t+1) = max(0, h_c(t) + dt/A*(q_in(t) - q_out_c(t)));
        theta_c(t) = 1;

    elseif (C_c(t) < climit)
        q_out_c(t) = co*Ao*sqrt(2*g*h_c(t))*min(1,h_c(t));
        h_c(t+1) = max(0, h_c(t) + dt/A*(q_in(t) - q_out_c(t)));
        theta_c(t) = 1;
    end
end
spill_c = 134285/10.764*max(0, h_c-10)/3.281/15/60.* C_c;

%% RBC-Both

h_b(1) = h0;
C_b(1) = C0;
% Operating policy
for t = 1:T
    % pond configuration
    elevation = [0, 2, 4, 6, 8, 10];
    area = [82971, 93258, 106100, 119152, 134285, 134285];
    A = interp1(elevation, area, h_b(t),'spline');

    C_b(t+1) = (C_b(t)*A*h_b(t)*exp(-k*dt)+c_in(t)*q_in(t)*dt)/(A*h_b(t)+q_in(t)*dt);

    if (C_b(t) > climit) && (h_b(t) < hlimit)
        q_out_b(t) = 0;
        h_b(t+1) = h_b(t) + dt/A*(q_in(t) - q_out_b(t));
        theta_b(t) = 0;

    elseif (h_b(t) >= hlimit)
        q_out_b(t) = q_desired;
        h_b(t+1) = max(0, h_b(t) + dt/A*(q_in(t) - q_out_b(t)));
        theta_b(t) = 1;

    elseif (C_b(t) < climit)
        if (h_b(t) < hlimit) && (h_b(t) >= hdes)
            q_out_b(t) = q_desired;
            h_b(t+1) = max(0, h_b(t) + dt/A*(q_in(t) - q_out_b(t)));
            theta_b(t) = q_desired/(co*Ao*sqrt(2*g*h_b(t)));
       
        elseif (h_b(t) < hdes)
            q_out_b(t) = co*Ao*sqrt(2*g*h_b(t))*min(1,h_b(t));
            h_b(t+1) = max(0, h_b(t) + dt/A*(q_in(t) - q_out_b(t)));
            theta_b(t) = 1;
        end

    end
end
spill_b = 134285/10.764*max(0, h_b-10)/3.281/15/60.* C_b;

%% Unit conversion from US to SI
qin_t = qin_t*0.028316846592;

nc.Xopt(:,1) = nc.Xopt(:,1)/3.281;
nc.Yopt(:,3) = nc.Yopt(:,3)*0.028316846592;

hlimit = hlimit/3.281;

falh = falh/3.281;
faly = faly*0.028316846592;

trueh = trueh/3.281;
truey = truey*0.028316846592;

mpc_ekfh = mpc_ekfh/3.281;
mpc_ekfy = mpc_ekfy*0.028316846592;

hh = hh/3.281;
yy = yy*0.028316846592;

q_out_q = q_out_q*0.028316846592;
q_out_c = q_out_c*0.028316846592;
q_out_b = q_out_b*0.028316846592;
h_q = h_q/3.281;
h_c = h_c/3.281;
h_b = h_b/3.281;

%% Perforamance comparison

% overflow
overflow = max(0,[(max(trueh) - 10/3.281)*134285/10.764, (max(h_c) - 10/3.281)*134285/10.764, (max(h_q) - 10/3.281)*134285/10.764, (max(h_b) - 10/3.281)*134285/10.764]);
overflow_percent = overflow/max(overflow)*100;

% peak outflow
outflow = [max(truey), max(q_out_c), max(q_out_q), max(q_out_b), max(nc.Yopt(:,3))];
outflow_percent = outflow/max(outflow)*100;

% cumulative load
load = [max(cumsum(truec.*truey)*10^(-3)*15*60), max(cumsum(C_c.*[0,q_out_c]+spill_c)*10^(-3)*15*60), max(cumsum(C_q.*[0,q_out_q]+spill_q)*10^(-3)*15*60), max(cumsum(C_b.*[0,q_out_b]+spill_b)*10^(-3)*15*60), max(cumsum(nc.Xopt(1:T,2).*nc.Yopt(1:T,3))*10^(-3)*15*60)];
load_percent = load/max(load)*100;

% control effort
coneff = [sum(diff(truemv).^2), sum(diff(theta_c).^2), sum(diff(theta_q).^2), sum(diff(theta_b).^2), sum(diff(nc.MVopt(1:T,1)).^2)];
coneff_percent = coneff/max(coneff)*100;

% outflow smoothness
smooth = [sum((truey - mean(truey)).^2), sum((q_out_c - mean(q_out_c)).^2), sum((q_out_q - mean(q_out_q)).^2), sum((q_out_b - mean(q_out_b)).^2), sum((nc.Yopt(1:T,3) - mean(nc.Yopt(1:T,3))).^2)];
smooth_percent = smooth/max(smooth)*100;

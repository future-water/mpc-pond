clc
clear
%% Initializing storm event runoff
load('Stormevent_Runoff.mat')

%%  Isolating the strom for analysis 

for i=1:800
    Temp1(i)=0;
end

Runoff=horzcat(transpose(A(1210:1737)),Temp1); % Extending the duration of dry period after storm in order to visualize the behavior of system after the event.

%% Generating Pollutant Loads 

for i=1:length(Runoff)
    if Runoff(i)==0
        Nitrate(i)=0;
    else 
        Nitrate(i)=( Runoff(i)*0.0004719474*100 + 3);
    end
end

%% Time Series for Simulink

Flow=timeseries(Runoff*0.0004719474*100); % Exporting runoff in cubic m per sec
Pollutant=timeseries(Nitrate); % Pollutant as mg per l

%% Plots 
figure('Units','centimeters',...
    'Position',[0 0 8.3 10],...
    'PaperPositionMode','auto');
subplot(2,1,1)
plot(Flow)
xlim([0 528])
title('Runoff')
ylabel('Discharge')

subplot(2,1,2)
plot(Pollutant)
title('Nitrate')
xlim([0 528])
ylabel('Concentration')

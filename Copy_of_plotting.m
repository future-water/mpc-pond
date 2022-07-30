fig = figure('units','normalized', 'outerposition', [0 0 1 1]);
c = [13 145 235]./255;
t = tiledlayout(4, 3, 'TileSpacing','compact', 'Padding','compact');
Nstep = length(orinc);

nexttile;
plot((1:Nstep), MD(1:Nstep,1), '-', 'LineWidth',2)
ylabel('Inflow (cfs)', 'FontSize',14,'FontWeight','bold')
xticklabels({})

nexttile(4);
plot((1:Nstep), MD(1:Nstep,2), '-','Color','r', 'LineWidth',2)
hold on
plot((1:Nstep), cin_t(1:Nstep,1), '-',	'Color','#EDB120','LineWidth',2)
hold off
legend({'True','Forecasted'})
ylabel('Inflow TSS (mg/L)', 'FontSize',14,'FontWeight','bold')
%xlim([0, 100])
xticklabels({})

nexttile;
nexttile(2,[2 1])
plot((1:Nstep)*Ts,nc.Xopt(1:Nstep,1),'--','Color','black','LineWidth',2)
hold on
plot((1:Nstep)*Ts, trueh, '-', 'LineWidth',2)
plot((1:Nstep)*Ts, orinh, '-', 'LineWidth',2)
plot((1:Nstep)*Ts, inith, '-', 'Color',c,'LineWidth',2)
hold off
legend({'No control','True','Forecasted','Control'})
ylabel('Pond height (ft)', 'FontSize',14,'FontWeight','bold')
ylim([0, 10])

nexttile;
nexttile(3,[2 1])
plot((1:Nstep)*Ts,nc.Yopt(1:Nstep,3),'--','Color','black','LineWidth',2)
hold on
plot((1:Nstep)*Ts, truey, '-', 'LineWidth',2)
plot((1:Nstep)*Ts, oriny, '-', 'LineWidth',2)
plot((1:Nstep)*Ts, inity, '-', 'Color',c,'LineWidth',2)
hold off
legend({'No control','True','Forecasted','Control'})
ylabel('Outflow (cfs)', 'FontSize',14,'FontWeight','bold')
xticklabels({})

nexttile;
nexttile(7,[2 1])
plot((1:Nstep)*Ts,nc.Xopt(1:Nstep,2),'--','Color','black','LineWidth',2)
hold on
plot((1:Nstep)*Ts, truec, '-','LineWidth',2)
plot((1:Nstep)*Ts, orinc, '-','LineWidth',2)
plot((1:Nstep)*Ts, initc, '-', 'Color',c,'LineWidth',2)
hold off
legend({'No control','True','Forecasted','Control'})
ylabel('TSS concentration (mg/L)','FontSize',14,'FontWeight','bold')

nexttile;
nexttile(8,[2 1])
plot((1:Nstep)*Ts,nc.MVopt(1:Nstep,1),'--','Color','black','LineWidth',2)
hold on
plot((1:Nstep)*Ts, truemv, '-', 'LineWidth',2)
plot((1:Nstep)*Ts, orinmv, '-', 'LineWidth',2)
plot((1:Nstep)*Ts, initmv, '-', 'Color',c,'LineWidth',2)
hold off
legend({'No control','True','Forecasted','Control'})
ylabel('Valve opening ratio (%)','FontSize',14,'FontWeight','bold')

nexttile;
nexttile(9,[2 1])
plot((1:Nstep)*Ts,cumsum(nc.Xopt(1:Nstep,2).*nc.Yopt(1:Nstep,3))*10^(-5)*15*60,'--','Color','black','LineWidth',2)
hold on
plot((1:Nstep)*Ts, cumsum(truec.*truey)*10^(-5)*15*60, '-', 'LineWidth',2)
plot((1:Nstep)*Ts, cumsum(orinc.*oriny)*10^(-5)*15*60, '-', 'LineWidth',2)
plot((1:Nstep)*Ts, cumsum(initc.*inity)*10^(-5)*15*60, '-', 'Color',c,'LineWidth',2)
hold off
legend({'No control','True','Forecasted','Control'})
ylabel('Cummulative Load (kg)','FontSize',14,'FontWeight','bold')

xlabel(t,'Time','FontSize',14,'FontWeight','bold')


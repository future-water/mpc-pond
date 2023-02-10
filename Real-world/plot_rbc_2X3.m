nIDs = 8;
alphabet = ('a':'z').';
chars = num2cell(alphabet(1:nIDs));
chars = chars.';
charlbl = strcat('(',chars,')'); % {'(a)','(b)','(c)','(d)'}

% precipitation data
dt = importdata('rainfall_15min_2021.csv');
prep = [zeros(48,1); dt.data(:,1)];

fig = figure('Units', 'points', 'position',[100 100 1500 900]);
c = [13 145 235]./255;
set(gcf(), 'DefaultAxesFontSize', 20)
tile = tiledlayout(9, 3, 'TileSpacing','compact', 'Padding','compact');
start = 24600;
start = 25220;
Nstep = 26275;


x = [0, start, start, 0];
y = [0, 0, 1.2*hlimit, 1.2*hlimit];

x2 = [Nstep, 34000, 34000, Nstep];
y2 = [0, 0, 1.2*hlimit, 1.2*hlimit];
ax = nexttile(1, [1,3])

plot((1:34000)*Ts,h_p(1:34000),'-','Color','black','LineWidth',2)
hold on
plot((1:34000)*Ts, hlimit*ones(1,34000), '--',  'color',[0.7 0.7 0.7], 'LineWidth',2)
plot((1:34000), h_c(1:34000),'-.','Color','#EDB120','LineWidth',2)
plot((1:34000), h_q(1:34000),':','Color','#D95319','LineWidth',2)
plot((1:34000), h_b(1:34000),'--','Color','#77AC30','LineWidth',2)
plot((1:34000)*Ts, mpc_ekfh(1:34000), '-', 'Color','#0072BD', 'LineWidth',3)
patch(x,y,[.7 .7 .7],'FaceAlpha',.5);
patch(x2,y2,[.7 .7 .7],'FaceAlpha',.5);
hold off
ylim([0, 1.2*hlimit])
xticks([20 2929 5617 8593 11473 14449 17329 20305 23281 26161 29137 32017])
xticklabels({'Jan.','Feb.','Mar.','Apr.','May','Jun.','Jul.','Aug.','Sep.','Oct.','Nov.','Dec.'})
xlim([1 34000])
text(0.005,0.8,charlbl{1},'Units','normalized','FontSize',24)

axh = nexttile(4, [2 1]);
plot((1:Nstep), qin_t(1:Nstep,1), '-','Color','#0000a7',  'LineWidth',3)
ylabel('Inflow (m^{3}/s)', 'FontSize',20,'FontWeight','bold')
ylim([0,2])
hold on
yyaxis right
axh.YDir = 'reverse';
axh.YColor = 'blue'; 
bar((start:Nstep), prep(start:Nstep,1).*25.4*4,'blue', 'BarWidth', 1)
ylabel('Precipitation (mm/hr)', 'FontSize',18,'FontWeight','bold', 'Color','blue')
ylim([0, 2*max(prep(start:Nstep,1).*25.4*4)])
xticklabels({})
xlim([start, Nstep])
text(0.020,0.9,charlbl{2},'Units','normalized','FontSize',24)

nexttile(10, [2 1]);
plot((1:Nstep), cin_t(1:Nstep,1), '-','Color','#c1272d','LineWidth',3)
hold on
plot((start:Nstep), MD(start:Nstep,2), ':','Color','#7E2F8E', 'LineWidth',3)
hold off
legend({'True pollutograph','False pollutograph'},'FontSize',16,'Box','off')
ylabel('TSS (mg/L)', 'FontSize',20,'FontWeight','bold')
xticklabels({})
xlim([start, Nstep])
ylim([0, 1.5*max(cin_t)])
text(0.020,0.9,charlbl{3},'Units','normalized','FontSize',24)


h1 = nexttile(5, [4,1])
plot((start:Nstep)*Ts,h_p(start:Nstep),'-','Color','black','LineWidth',4)
hold on
plot((start:Nstep)*Ts, hlimit*ones(1,Nstep - start+1), '--',  'color',[0.7 0.7 0.7], 'LineWidth',4)
plot((start:Nstep), h_c(start:Nstep),'-.','Color','#EDB120','LineWidth',4)
plot((start:Nstep), h_q(start:Nstep),':','Color','#D95319','LineWidth',4)
plot((start:Nstep), h_b(start:Nstep),'--','Color','#77AC30','LineWidth',4)
plot((start:Nstep)*Ts, mpc_ekfh(start:Nstep), '-', 'Color','#0072BD', 'LineWidth',5)
hold off
legend({'','Maximum allowable height'},'FontSize',16,'Box','off')
ylabel('Pond height (m)', 'FontSize',20,'FontWeight','bold')
ylim([0, 1.25*hlimit])
xticklabels({})
xlim([start, Nstep])
text(0.020,0.95,charlbl{4},'Units','normalized','FontSize',24)


nexttile;
nexttile(6, [4 1])
plot((start:Nstep), q_out_c_cut(start:Nstep),'-.','Color','#EDB120','LineWidth',3)
hold on
plot((start:Nstep)*Ts,q_out_p(start:Nstep),'-','Color','black','LineWidth',4)
plot((start:Nstep), q_out_b(start:Nstep),'--','Color','#77AC30','LineWidth',3)
plot((start:Nstep), q_out_q(start:Nstep),':','Color','#D95319','LineWidth',4)
plot((start:Nstep)*Ts, mpc_ekfy(start:Nstep), '-', 'Color','#0072BD', 'LineWidth',5)
hold off
ylabel('Outflow (m^{3}/s)', 'FontSize',20,'FontWeight','bold')
xticklabels({})
xlim([start, Nstep])
text(0.020,0.95,charlbl{5},'Units','normalized','FontSize',24)


nexttile;
nexttile(16,[4 1])
plot((start:Nstep)*Ts,C_p(start:Nstep),'-','Color','black','LineWidth',4)
hold on
plot((start:Nstep), C_c(start:Nstep),'-.','Color','#EDB120','LineWidth',4)
plot((start:Nstep), C_q(start:Nstep),':','Color','#D95319','LineWidth',4)
plot((start:Nstep), C_b(start:Nstep),'--','Color','#77AC30','LineWidth',4)
plot((start:Nstep)*Ts, mpc_ekfc(start:Nstep), '-', 'Color','#0072BD','LineWidth',5)
hold off
ylabel('TSS concentration (mg/L)','FontSize',20,'FontWeight','bold')
xlim([start, Nstep])
xticks([24770 25250 25730 26210])
xticklabels({'Sep. 17','Sep. 22','Sep. 27','Oct. 1'})
ylim([0, 1.2*max(mpc_ekfc(start:Nstep))])
text(0.020,0.95,charlbl{6},'Units','normalized','FontSize',24)


nexttile;
nexttile(17,[4 1])
plot((start:Nstep)*Ts,ones(Nstep - start +1 ,1),'-','Color','black','LineWidth',4)
hold on
plot((start:Nstep), theta_c(start:Nstep),'-.','Color','#EDB120','LineWidth',3)
plot((start:Nstep), theta_b(start:Nstep),'--','Color','#77AC30','LineWidth',3)
plot((start:Nstep), theta_q(start:Nstep),':','Color','#D95319','LineWidth',4)
plot((start:Nstep)*Ts, mpc_ekfmv(start:Nstep), '-', 'Color','#0072BD', 'LineWidth',5)
hold off
ylabel('Valve opening ratio','FontSize',20,'FontWeight','bold')
xlim([start, Nstep])
xticks([24770 25250 25730 26210])
xticklabels({'Sep. 17','Sep. 22','Sep. 27','Oct. 1'})
ylim([0,1])
text(0.020,0.95,charlbl{7},'Units','normalized','FontSize',24)


nexttile;
nexttile(18,[4 1])
plot((start:Nstep)*Ts,cumsum(C_p(start:Nstep).*q_out_p(start:Nstep))*10^(-3)*15*60,'-','Color','black','LineWidth',4)
hold on
plot((start:Nstep)*Ts, cumsum(mpc_ekfc(start:Nstep).*mpc_ekfy(start:Nstep))*10^(-3)*15*60, '-', 'Color','#0072BD', 'LineWidth',5)
plot((start:Nstep), cumsum(C_c(start:Nstep).*q_out_c_cut(start:Nstep)+spill_c(start:Nstep))*10^(-3)*15*60,'-.','Color','#EDB120','LineWidth',4)
plot((start:Nstep), cumsum(C_q(start:Nstep).*q_out_q(start:Nstep)+spill_q(start:Nstep))*10^(-3)*15*60,':','Color','#D95319','LineWidth',4)
plot((start:Nstep), cumsum(C_b(start:Nstep).*q_out_b(start:Nstep)+spill_b(start:Nstep))*10^(-3)*15*60,'--','Color','#77AC30','LineWidth',4)
hold off
ylabel('Cummulative Load (kg)','FontSize',20,'FontWeight','bold')
xlim([start, Nstep])
xticks([24770 25250 25730 26210])
xticklabels({'Sep. 17','Sep. 22','Sep. 27','Oct. 1'})
text(0.020,0.95,charlbl{8},'Units','normalized','FontSize',24)
leg = legend({'Passive','MPC-EKF','RBC-Conc','RBC-Outflow','RBC-Both'},'Location','southoutside','NumColumns',5, 'FontSize',24,'FontWeight','bold');
leg.Layout.Tile = 'north';

MagInset(tile, h1, [25300 25550 2.85 3.45], [25900 26250 1.7 2.9], {'NE','NW';'SE','SW'});
xticklabels({})
nIDs = 7;
alphabet = ('a':'z').';
chars = num2cell(alphabet(1:nIDs));
chars = chars.';
charlbl = strcat('(',chars,')'); % {'(a)','(b)','(c)','(d)'}

% precipitation data
dt = importdata('rainfall_15min_2021.csv');
prep = dt.data(:,1);

c = [13 145 235]./255;
set(gcf(), 'DefaultAxesFontSize', 16)
set(gcf, 'Position', get(0, 'Screensize'));
tile = tiledlayout(12, 1, 'TileSpacing','compact', 'Padding','compact');
start = 1;
Nstep = 34000;

axh = nexttile(1);
plot((1:Nstep), qin_t(1:Nstep,1), '-','Color','#0000a7',  'LineWidth',3)
ylabel('Inflow (m^{3}/s)', 'FontSize',14,'FontWeight','bold')
ylim([0,5])
hold on
yyaxis right
axh.YDir = 'reverse';
axh.YColor = 'blue'; 
bar((start:Nstep), prep(start:Nstep,1).*25.4*4,'blue', 'BarWidth', 10)
ylabel('P (mm/hr)', 'FontSize',14,'FontWeight','bold', 'Color','blue')
xticklabels({})
xlim([start, Nstep])
text(0.020,0.8,charlbl{1},'Units','normalized','FontSize',20)

nexttile(2);
plot((1:Nstep), cin_t(1:Nstep,1), '-','Color','#c1272d','LineWidth',3)
hold on
plot((start:Nstep), MD(start:Nstep,2), ':','Color','#7E2F8E', 'LineWidth',3)
hold off
legend({'True pollutograph','False pollutograph'},'FontSize',14,'Box','off','NumColumns',2)
ylabel('TSS (mg/L)', 'FontSize',14,'FontWeight','bold')
xticklabels({})
xlim([start, Nstep])
ylim([0, 1.5*max(cin_t)])
text(0.020,0.8,charlbl{2},'Units','normalized','FontSize',20)


nexttile(3, [2,1])
plot((start:Nstep)*Ts,h_p(start:Nstep),'-','Color','black','LineWidth',2)
hold on
plot((start:Nstep)*Ts, hlimit*ones(1,Nstep - start+1), '--',  'color',[0.7 0.7 0.7], 'LineWidth',2)
plot((start:Nstep), h_c(start:Nstep),'-.','Color','#EDB120','LineWidth',2)
plot((start:Nstep), h_q(start:Nstep),':','Color','#D95319','LineWidth',2)
plot((start:Nstep), h_b(start:Nstep),'--','Color','#77AC30','LineWidth',2)
plot((start:Nstep)*Ts, mpc_ekfh(start:Nstep), '-', 'Color','#0072BD', 'LineWidth',3)
hold off
legend({'','Maximum allowable height'},'FontSize',14,'Box','off')
ylabel('Pond height (m)', 'FontSize',14,'FontWeight','bold')
ylim([0, 1.2*hlimit])
xticklabels({})
xlim([start, Nstep])
text(0.020,0.9,charlbl{3},'Units','normalized','FontSize',20)


nexttile;
nexttile(5, [2 1])
plot((start:Nstep)*Ts,q_out_p(start:Nstep),'-','Color','black','LineWidth',2)
hold on
plot((start:Nstep), q_out_c_cut(start:Nstep),'-.','Color','#EDB120','LineWidth',2)
plot((start:Nstep), q_out_q(start:Nstep),':','Color','#D95319','LineWidth',2)
plot((start:Nstep), q_out_b(start:Nstep),'--','Color','#77AC30','LineWidth',2)
plot((start:Nstep)*Ts, mpc_ekfy(start:Nstep), '-', 'Color','#0072BD', 'LineWidth',3)
hold off
ylabel('Outflow (m^{3}/s)', 'FontSize',14,'FontWeight','bold')
xticklabels({})
xlim([start, Nstep])
text(0.020,0.9,charlbl{4},'Units','normalized','FontSize',20)


nexttile;
nexttile(7,[2 1])
plot((start:Nstep)*Ts,C_p(start:Nstep),'-','Color','black','LineWidth',2)
hold on
plot((start:Nstep), C_c(start:Nstep),'-.','Color','#EDB120','LineWidth',2)
plot((start:Nstep), C_q(start:Nstep),':','Color','#D95319','LineWidth',2)
plot((start:Nstep), C_b(start:Nstep),'--','Color','#77AC30','LineWidth',2)
plot((start:Nstep)*Ts, mpc_ekfc(start:Nstep), '-', 'Color','#0072BD','LineWidth',3)
hold off
ylabel('TSS conc. (mg/L)','FontSize',14,'FontWeight','bold')
xlim([start, Nstep])
xticklabels({})
ylim([0, 1.2*max(mpc_ekfc(start:Nstep))])
text(0.020,0.9,charlbl{5},'Units','normalized','FontSize',20)


nexttile;
nexttile(9,[2 1])
plot((start:Nstep)*Ts,ones(Nstep - start +1 ,1),'-','Color','black','LineWidth',2)
hold on
plot((start:Nstep), theta_c(start:Nstep),'-.','Color','#EDB120','LineWidth',2)
plot((start:Nstep), theta_q(start:Nstep),':','Color','#D95319','LineWidth',2)
plot((start:Nstep), theta_b(start:Nstep),'--','Color','#77AC30','LineWidth',2)
plot((start:Nstep)*Ts, mpc_ekfmv(start:Nstep), '-', 'Color','#0072BD', 'LineWidth',3)
hold off
ylabel('Valve opening ratio','FontSize',14,'FontWeight','bold')
xticklabels({})
xlim([start, Nstep])
ylim([0,1])
text(0.020,0.9,charlbl{6},'Units','normalized','FontSize',20)


nexttile;
nexttile(11,[2 1])
plot((start:Nstep)*Ts,cumsum(C_p(start:Nstep).*q_out_p(start:Nstep))*10^(-3)*15*60,'-','Color','black','LineWidth',2)
hold on
plot((start:Nstep)*Ts, cumsum(mpc_ekfc(start:Nstep).*mpc_ekfy(start:Nstep))*10^(-3)*15*60, '-', 'Color','#0072BD', 'LineWidth',3)
plot((start:Nstep), cumsum(C_c(start:Nstep).*q_out_c_cut(start:Nstep)+spill_c(start:Nstep))*10^(-3)*15*60,'-.','Color','#EDB120','LineWidth',2)
plot((start:Nstep), cumsum(C_q(start:Nstep).*q_out_q(start:Nstep)+spill_q(start:Nstep))*10^(-3)*15*60,':','Color','#D95319','LineWidth',2)
plot((start:Nstep), cumsum(C_b(start:Nstep).*q_out_b(start:Nstep)+spill_b(start:Nstep))*10^(-3)*15*60,'--','Color','#77AC30','LineWidth',2)
hold off
ylabel('Cummulative Load (kg)','FontSize',14,'FontWeight','bold')
xlim([start, Nstep])
xticks([0 2929 5617 8593 11473 14449 17329 20305 23281 26161 29137 32017])
xticklabels({'Jan.','Feb.','Mar.','Apr.','May','Jun.','Jul.','Aug.','Sep.','Oct.','Nov.','Dec.'})
text(0.020,0.9,charlbl{7},'Units','normalized','FontSize',20)
leg = legend({'Passive','MPC-EKF','RBC-Conc','RBC-Outflow','RBC-Both'},'Location','southoutside','NumColumns',5, 'FontSize',14,'FontWeight','bold');
leg.Layout.Tile = 'north';
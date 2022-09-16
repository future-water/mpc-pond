nIDs = 7;
alphabet = ('a':'z').';
chars = num2cell(alphabet(1:nIDs));
chars = chars.';
charlbl = strcat('(',chars,')'); % {'(a)','(b)','(c)','(d)'}

fig = figure('Units', 'points', 'position',[100 100 1500 900]);
c = [13 145 235]./255;

tile = tiledlayout(4, 3, 'TileSpacing','compact', 'Padding','compact');
Nstep = 820;

nexttile;
plot((1:Nstep), qin_t(1:Nstep,1), '-','Color','#0000a7',  'LineWidth',3)
ylabel('Inflow (m^{3}/s)', 'FontSize',14,'FontWeight','bold')
xticklabels({})
xlim([0, Nstep])
text(0.010,0.90,charlbl{1},'Units','normalized','FontSize',14)

nexttile(4);
plot((1:Nstep), cin_t(1:Nstep,1), '-','Color','#c1272d','LineWidth',3)
hold on
plot((1:Nstep), MD(1:Nstep,2), ':','Color','#77AC30', 'LineWidth',3)
hold off
legend({'True pollutograph','False pollutograph'},'FontSize',12,'Box','off')
ylabel('Inflow TSS (mg/L)', 'FontSize',14,'FontWeight','bold')
xticklabels({})
xlim([0, Nstep])
ylim([0, 1.3*max(cin_t)])
text(0.010,0.90,charlbl{2},'Units','normalized','FontSize',14)


nexttile(2,[2 1])
plot((1:Nstep)*Ts,nc.Xopt(1:Nstep,1),'-','Color','black','LineWidth',2)
hold on
plot((1:Nstep)*Ts, trueh(1:Nstep), '-', 'Color','#0072BD', 'LineWidth',4)
plot((1:Nstep)*Ts, hlimit*ones(1,Nstep), '--',  'color',[0.7 0.7 0.7], 'LineWidth',2)
plot(h_c,'-.','Color','#EDB120','LineWidth',3)
plot(h_q,':','Color','#D95319','LineWidth',3)
hold off
ylabel('Pond height (m)', 'FontSize',14,'FontWeight','bold')
ylim([0, 1.2*hlimit])
xticklabels({})
xlim([0, Nstep])
text(0.010,0.95,charlbl{3},'Units','normalized','FontSize',14)


nexttile;
nexttile(3, [2 1])
plot((1:Nstep)*Ts,nc.Yopt(1:Nstep,3),'-','Color','black','LineWidth',2)
hold on
plot((1:Nstep)*Ts, truey(1:Nstep), '-', 'Color','#0072BD', 'LineWidth',4)
plot(q_out_c(1:Nstep),'-.','Color','#EDB120','LineWidth',3)
plot(q_out_q(1:Nstep),':','Color','#D95319','LineWidth',3)
hold off
ylabel('Outflow (m^{3}/s)', 'FontSize',14,'FontWeight','bold')
xticklabels({})
xlim([0, Nstep])
text(0.010,0.95,charlbl{4},'Units','normalized','FontSize',14)


nexttile;
nexttile(7,[2 1])
plot((1:Nstep)*Ts,nc.Xopt(1:Nstep,2),'-','Color','black','LineWidth',2)
hold on
plot((1:Nstep)*Ts, truec(1:Nstep), '-', 'Color','#0072BD','LineWidth',4)
plot(C_c(1:Nstep),'-.','Color','#EDB120','LineWidth',3)
plot(C_q(1:Nstep),':','Color','#D95319','LineWidth',3)
hold off
ylabel('TSS concentration (mg/L)','FontSize',14,'FontWeight','bold')
xticks([0 24 48 72 96 120 144 168 192]*4)
xticklabels({'0','24','48','72','96','120','144','168','192'})
xlim([0, Nstep])
ylim([0 27])
text(0.010,0.95,charlbl{5},'Units','normalized','FontSize',14)


nexttile;
nexttile(8,[2 1])
plot((1:Nstep)*Ts,nc.MVopt(1:Nstep,1),'-','Color','black','LineWidth',2)
hold on
plot((1:Nstep)*Ts, truemv(1:Nstep), '-', 'Color','#0072BD', 'LineWidth',4)
plot(theta_c(1:Nstep),'-.','Color','#EDB120','LineWidth',3)
plot(theta_q(1:Nstep),':','Color','#D95319','LineWidth',3)
hold off
ylabel('Valve opening ratio (%)','FontSize',14,'FontWeight','bold')
xlim([0, Nstep])
xticks([0 24 48 72 96 120 144 168 192]*4)
xticklabels({'0','24','48','72','96','120','144','168','192'})
ylim([0,1])
text(0.010,0.95,charlbl{6},'Units','normalized','FontSize',14)
leg = legend({'Passive','MPC','RBC-Conc.','RBC-Outflow'},'Location','southoutside','NumColumns',5, 'FontSize',14','FontWeight','bold');
leg.Layout.Tile = 'north';

nexttile;
h2 = nexttile(9,[2 1])
plot((1:Nstep)*Ts,cumsum(nc.Xopt(1:Nstep,2).*nc.Yopt(1:Nstep,3))*10^(-3)*15*60,'-','Color','black','LineWidth',2)
hold on
plot((1:Nstep)*Ts, cumsum(truec(1:Nstep).*truey(1:Nstep))*10^(-3)*15*60, '-', 'Color','#0072BD', 'LineWidth',4)
plot(cumsum(C_c(1:Nstep).*q_out_c(1:Nstep)+spill_c(1:Nstep))*10^(-3)*15*60,'-.','Color','#EDB120','LineWidth',3)
plot(cumsum(C_q(1:Nstep).*q_out_q(1:Nstep)+spill_q(1:Nstep))*10^(-3)*15*60,':','Color','#D95319','LineWidth',3)
hold off
ylabel('Cummulative Load (kg)','FontSize',14,'FontWeight','bold')
xlabel(tile, 'Time (h)','FontSize',14,'FontWeight','bold')
xlim([0, Nstep])
xticks([0 24 48 72 96 120 144 168 192]*4)
xticklabels({'0','24','48','72','96','120','144','168','192'})
text(0.010,0.95,charlbl{7},'Units','normalized','FontSize',14)


MagInset(tile, h2, [700 750 80 800], [480 780 3000 7000], {'NW','SW';'NE','SE'});
xticks([0 24 48 72 96 120 144 168 192]*4)
xticklabels({'0','24','48','72','96','120','144','168','192'})

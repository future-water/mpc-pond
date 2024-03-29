nIDs = 7;
alphabet = ('a':'z').';
chars = num2cell(alphabet(1:nIDs));
chars = chars.';
charlbl = strcat('(',chars,')'); % {'(a)','(b)','(c)','(d)'}

fig = figure('Units', 'points', 'position',[100 100 750 1000]);
c = [13 145 235]./255;

tile = tiledlayout(6, 2, 'TileSpacing','compact', 'Padding','compact');
Nstep = 820;

nexttile;
plot((1:Nstep), qin_t(1:Nstep,1), '-','Color','#0000a7',  'LineWidth',3)
ylabel('Inflow (m^{3}/s)', 'FontSize',14,'FontWeight','bold')
xticklabels({})
xlim([0, Nstep])
text(0.010,0.90,charlbl{1},'Units','normalized','FontSize',14)

nexttile(3);
plot((1:Nstep), cin_t(1:Nstep,1), '-','Color','#c1272d','LineWidth',3)
hold on
plot((1:Nstep), cin_f(1:Nstep,:), ':', 'LineWidth',3)
hold off
legend({'True pollutograph','False pollutograph'},'FontSize',12,'Box','off')
ylabel('Inflow TSS (mg/L)', 'FontSize',14,'FontWeight','bold')
xticklabels({})
xlim([0, Nstep])
ylim([0 45])
text(0.010,0.90,charlbl{2},'Units','normalized','FontSize',14)

nexttile(2,[2 1])
plot((1:Nstep)*Ts,nc.Xopt(1:Nstep,1),'-','Color','black','LineWidth',2)
hold on
plot((1:Nstep)*Ts, trueh(1:Nstep), '-', 'Color','#0072BD', 'LineWidth',2)
plot((1:Nstep)*Ts, falhh(1:Nstep,:), ':','Color','#77AC30','LineWidth',2)
plot((1:Nstep)*Ts, hh(1:Nstep,:), '-','LineWidth',3)
plot((1:Nstep)*Ts, hlimit*ones(1,Nstep), '--',  'color',[0.7 0.7 0.7], 'LineWidth',2)
hold off
ylabel('Pond height (m)', 'FontSize',14,'FontWeight','bold')
ylim([0, 1.1*hlimit])
xticklabels({})
xlim([0, Nstep])
text(0.010,0.95,charlbl{3},'Units','normalized','FontSize',14)

nexttile;
nexttile(5, [2 1])
plot((1:Nstep)*Ts,nc.Yopt(1:Nstep,3),'-','Color','black','LineWidth',2)
hold on
plot((1:Nstep)*Ts, truey(1:Nstep), '-', 'Color','#0072BD', 'LineWidth',2)
plot((1:Nstep)*Ts, falyy(1:Nstep,:), ':', 'Color','#77AC30','LineWidth',2)
plot((1:Nstep)*Ts, yy(1:Nstep,:), '-', 'Color','#4DBEEE','LineWidth',3)
hold off
ylabel('Outflow (m^{3}/s)', 'FontSize',14,'FontWeight','bold')
xticklabels({})
xlim([0, Nstep])
text(0.010,0.95,charlbl{4},'Units','normalized','FontSize',14)

nexttile;
h1 = nexttile(6,[2 1])
plot((1:Nstep)*Ts,nc.Xopt(1:Nstep,2),'-','Color','black','LineWidth',2)
hold on
plot((1:Nstep)*Ts, truec(1:Nstep), '-', 'Color','#0072BD','LineWidth',2)
plot((1:Nstep)*Ts, falcc(1:Nstep,:), ':','Color','#77AC30','LineWidth',2)
plot((1:Nstep)*Ts, cc(1:Nstep,:), '-', 'Color','#4DBEEE','LineWidth',3)
hold off
ylabel('TSS concentration (mg/L)','FontSize',14,'FontWeight','bold')
xticklabels({})
xlim([0, Nstep])
ylim([0, 37])
text(0.010,0.95,charlbl{5},'Units','normalized','FontSize',14)

nexttile;
h3 = nexttile(9,[2 1])
plot((1:Nstep)*Ts,nc.MVopt(1:Nstep,1),'-','Color','black','LineWidth',2)
hold on
plot((1:Nstep)*Ts, truemv(1:Nstep), '-', 'Color','#0072BD', 'LineWidth',2)
plot((1:Nstep)*Ts, falmvmv(1:Nstep,:), ':','Color','#77AC30', 'LineWidth',2)
plot((1:Nstep)*Ts, mvmv(1:Nstep,:), '-', 'Color',c,'LineWidth',3)
hold off
ylabel('Valve opening ratio (%)','FontSize',14,'FontWeight','bold')
xlim([0, Nstep])
xticks([0 24 48 72 96 120 144 168 192]*4)
xticklabels({'0','24','48','72','96','120','144','168','192'})
ylim([0,1])
text(0.010,0.95,charlbl{6},'Units','normalized','FontSize',14)
leg = legend({'Passive','MPC-EKF with noise','','','','','','','','','','','','','','','','','','','','MPC-True','MPC-False','MPC-EKF'},'Location','southoutside','NumColumns',5, 'FontSize',14','FontWeight','bold');
leg.Layout.Tile = 'north';

nexttile;
h2 = nexttile(10,[2 1])
plot((1:Nstep)*Ts,cumsum(nc.Xopt(1:Nstep,2).*nc.Yopt(1:Nstep,3))*10^(-3)*15*60,'-','Color','black','LineWidth',2)
hold on
text(0.010,0.95,charlbl{7},'Units','normalized','FontSize',14)
plot((1:Nstep)*Ts, cumsum(truec(1:Nstep).*truey(1:Nstep))*10^(-3)*15*60, '-', 'Color','#0072BD', 'LineWidth',2)
plot((1:Nstep)*Ts, cumsum(truec(1:Nstep,:).*falyy(1:Nstep,:))*10^(-3)*15*60, ':','Color','#77AC30', 'LineWidth',2)
plot((1:Nstep)*Ts, cumsum(cc(1:Nstep,:).*yy(1:Nstep,:))*10^(-3)*15*60, '-', 'Color','#4DBEEE','LineWidth',3)
hold off
ylabel('Cummulative Load (kg)','FontSize',14,'FontWeight','bold')
xlabel(tile, 'Time (h)','FontSize',14,'FontWeight','bold')
xlim([0, Nstep])

xticks([0 24 48 72 96 120 144 168 192]*4)
xticklabels({'0','24','48','72','96','120','144','168','192'})

MagInset(tile, h1, [10 100 22 33], [400 750 10 30], {'NE','NW';'SE','SW'});
xticks([0 12 24 48 72 96 120 144 168 192]*4)
xticklabels({'0','12','24','48','72','96','120','144','168','192'})

MagInset(tile, h2, [750 800 40 80], [450 700 150 400], {'NW','SW';'NE','SE'});
xticks([0 24 48 72 96 120 144 168 192]*4)
xticklabels({'0','24','48','72','96','120','144','168','192'})

MagInset(tile, h3, [0 100 0 0.02], [150 350 0.7 0.95], {'NW','SW';'NE','SE'});
xticks([0 12 24 48 72 96 120 144 168 192]*4)
xticklabels({'0','12','24','48','72','96','120','144','168','192'})


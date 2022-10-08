nIDs = 7;
alphabet = ('a':'z').';
chars = num2cell(alphabet(1:nIDs));
chars = chars.';
charlbl = strcat('(',chars,')'); % {'(a)','(b)','(c)','(d)'}

ctable1 = [jet(8), 0.5*ones(8,1)];
ctable2 = jet(8);

fig = figure('Units', 'points', 'position',[100 100 1500 450]);
c = [13 145 235]./255;

tile = tiledlayout(1, 25, 'TileSpacing','compact', 'Padding','compact');
Nstep = 820;

nexttile;
nexttile(1, [1 8])
plot((1:Nstep)*Ts,nc.Yopt(1:Nstep,3),'-','Color','black','LineWidth',2)
hold on
for i = 1:8
    plot((1:Nstep)*Ts, falyy(1:Nstep,i), ':','Color',ctable1(i,:), 'LineWidth',4)
end
for i = 1:8
    plot((1:Nstep)*Ts, yy(1:Nstep,i), '-','Color',ctable2(i,:), 'LineWidth',2)
end
plot((1:Nstep)*Ts, truey(1:Nstep), '-', 'Color','#0072BD', 'LineWidth',4)
hold off
ylabel('Outflow (m^{3}/s)', 'FontSize',14,'FontWeight','bold')
xticks([0 24 48 72 96 120 144 168 192]*4)
xticklabels({'0','24','48','72','96','120','144','168','192'})
xlim([0, Nstep])
text(0.010,0.95,charlbl{1},'Units','normalized','FontSize',14)

nexttile;
h3 = nexttile(9, [1 8])
plot((1:Nstep)*Ts,nc.MVopt(1:Nstep,1),'-','Color','black','LineWidth',2, 'DisplayName','cos(x)')
hold on
for i = 1:8
    plot((1:Nstep)*Ts, falmvmv(1:Nstep,i), ':','Color',ctable1(i,:), 'LineWidth',4)
end
for i = 1:8
    plot((1:Nstep)*Ts, mvmv(1:Nstep,i), '-','Color',ctable2(i,:), 'LineWidth',2)
end
plot((1:Nstep)*Ts, truemv(1:Nstep), '-.', 'Color',[1 1 1 0.8], 'LineWidth',3)
hold off
ylabel('Valve opening ratio (%)','FontSize',14,'FontWeight','bold')
xlim([0, Nstep])
xticks([0 24 48 72 96 120 144 168 192]*4)
xticklabels({'0','24','48','72','96','120','144','168','192'})
ylim([0,1])
text(0.010,0.95,charlbl{2},'Units','normalized','FontSize',14)
leg = legend({'Passive','MPC','','','','','','','','','MPC-EKF','','','','','','','','',''},'Location','southoutside','NumColumns',5, 'FontSize',14','FontWeight','bold');
leg.Layout.Tile = 'north';

nexttile;
h2 = nexttile(17, [1 8])
plot((1:Nstep)*Ts,cumsum(nc.Xopt(1:Nstep,2).*nc.Yopt(1:Nstep,3))*10^(-3)*15*60,'-','Color','black','LineWidth',2)
hold on
text(0.010,0.95,charlbl{3},'Units','normalized','FontSize',14)
for i = 1:8
    plot((1:Nstep)*Ts, cumsum(faltruecc(1:Nstep,i).*falyy(1:Nstep,i))*10^(-3)*15*60, ':','Color',ctable1(i,:), 'LineWidth',4)
end
for i = 1:8
    plot((1:Nstep)*Ts, cumsum(cc(1:Nstep,i).*yy(1:Nstep,i))*10^(-3)*15*60, '-','Color',ctable2(i,:), 'LineWidth',4)
end
plot((1:Nstep)*Ts, cumsum(truec(1:Nstep).*truey(1:Nstep))*10^(-3)*15*60, '-', 'Color','#0072BD', 'LineWidth',2)
hold off
ylabel('Cummulative Load (kg)','FontSize',14,'FontWeight','bold')
xlabel(tile, 'Time (h)','FontSize',14,'FontWeight','bold')
xlim([0, Nstep])

xticks([0 24 48 72 96 120 144 168 192]*4)
xticklabels({'0','24','48','72','96','120','144','168','192'})

colormap(jet(8))
cb = colorbar('FontSize',12,'fontweight','bold','position',[0.6767+.255 0.0891 .02 0.82], 'Location','eastoutside', 'Ticks',linspace(0,1,17),...
         'TickLabels',{'','25%','','50%','','75%','','100%','','125%','','150%','','175%','','200%',''})
title(cb, 'EMC scenarios');
cb.Box = 'off';
cb.TickLength = 0;

MagInset(tile, h2, [750 800 55 65], [450 700 250 400], {'NW','SW';'NE','SE'});
xticks([0 24 48 72 96 120 144 168 192]*4)
xticklabels({'0','24','48','72','96','120','144','168','192'})

MagInset(tile, h3, [0 100 0 0.07], [550 800 0.1 0.35], {'NW','NW';'NE','SE'});
xticks([0 12 24 48 72 96 120 144 168 192]*4)
xticklabels({'0','12','24','48','72','96','120','144','168','192'})
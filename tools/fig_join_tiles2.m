
clear;
close all;

load('color_purple_to_green.mat')
leg = {'-2 %/m','-1.5 %/m','-1 %/m','-0.5 %/m','0 %/m','0.5 %/m','1 %/m','1.5 %/m','2 %/m'};


fig1 = openfig(['gradsim_paper/energy_gain/fig/energygain.fig']);
oax1 = fig1.Children(2);

fig2 = openfig(['gradsim_paper/amplitude/fig/maxlongxi6.fig']);
oax2 = fig2.Children(2);



figx = figure(10);
figx.Units = 'pixels';
figx.OuterPosition = [100 100 800 400];
colororder(cc);

tt = tiledlayout(1,2);
tt.TileSpacing = 'compact';
tt.Padding = 'compact';

ax(1) = nexttile;

copyobj(oax1.Children,ax(1));
xlim([0 10])
ylim([-50 2000])
ylabel('energy gain (MeV)')
ll = legend(ax(1),leg,'Location','northwest','AutoUpdate','off');
legend('boxoff')

% xline(4,'--','LineWidth',1,'color',[0 0.4470 0.7410]);

ax(2) = nexttile;
colororder(ax(2),cc(5:end,:));

copyobj(oax2.Children,ax(2));
% xline(4,'--','LineWidth',1,'color',[0 0.4470 0.7410]);
xlim([0 10.0])
ylabel('maximum E_z (MV/m)')

xlabel(tt,'z (m)')

P = Plotty();
P.plots_dir = 'gradsim_paper/join/';
P.plot_name = 'energygainamplitude';
P.fig_handle = figx;
P.save_plot();


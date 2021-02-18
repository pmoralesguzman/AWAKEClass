
clear;
close all;

load('color_red_to_blue.mat');

fig1 = openfig(['gradsim_paper/dephasing_lines/fig/dephasing1471allgrads0.fig']);
oax1 = fig1.Children;

fig2 = openfig(['gradsim_paper/chargeinxi/fig/chargeinxi.fig']);
oax2 = fig2.Children;



figx = figure(10);
figx.Units = 'pixels';
figx.OuterPosition = [100 100 1200 800];
colororder(ccrb);

tt = tiledlayout(2,3);
tt.TileSpacing = 'compact';
tt.Padding = 'compact';

ngrid = fig1.Children.GridSize(2);

for t = 1:fig1.Children.GridSize(2)

ax(t) = nexttile(3-t+1);

copyobj(oax1.Children(1+t).Children(:),ax(t));
copyobj(oax1.Children(1+t).Title,ax(t));


% xline(4,'--','LineWidth',1,'color',[0 0.4470 0.7410]);

ax(ngrid+t) = nexttile(2*ngrid-t+1);
% colororder(ax(2),cc(5:end,:));

copyobj(oax2.Children(1+t).Children(:),ax(ngrid+t));
% xline(4,'--','LineWidth',1,'color',[0 0.4470 0.7410]);
% xlim([0 10.0])
% ylabel('maximum E_z (MV/m)')

% xlabel(tt,'z (m)')

end


xlabel(tt,'z (m)')

P = Plotty();
P.plots_dir = 'gradsim_paper/join/';
P.plot_name = 'dephasing_charge';
P.fig_handle = figx;
P.save_plot();


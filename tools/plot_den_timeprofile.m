%________________________________________________________________________
% Script to manually (not using Plotty) plot the density profiles of the
% simulations (time profiles)
%
% Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 07/12/2020
%________________________________________________________________________


% close all;
clear;

% file location variables
datadirs_sim = {'g0'};

dataformat = 'mat';
useAvg = false;
dump = 100;

% saving data
save_flag = true;
save_format = {'pdf','png'};

% plasma properties
plasmaden = 1.81e14;

% choose fields to plot
wakefields_direction = ''; % trans, long

% choose species density to plot
species = 'proton_beam';

% choose limits (in cm for transverse, in ps for longitudinal)
trans_range = [0 0.16];
xi_range = [-1,618]; %ps

% choose property to plot
property_plot = 'density'; % density, wakefields, both

% create movie or not
create_movie = false;

% choose between normalized and denormalized units
denormalize_flag = true; % true, false

% choose scaling
plot_scale = 'linear';

% choose if make pause or not
make_pause = false;

% figure number
fig_number = 1;

letters = {'g = 0 %/m'};


% directory to save the plots
plots_dir = ['gradsim/','field_density_longprofile/'];

P = Plotty('datadir',datadirs_sim{1},'dataformat',dataformat,...
    'useAvg',useAvg,...
    'plot_scale',plot_scale,...
    'save_flag',save_flag,'save_format',save_format,'plots_dir',plots_dir,...
    'create_movie',create_movie,...
    'plasmaden',plasmaden,'trans_range',trans_range,'xi_range',xi_range,...
    'wakefields_direction',wakefields_direction,'species',species,...
    'property_plot',property_plot,'denormalize_flag',denormalize_flag,...
    'make_pause',make_pause,'fig_number',fig_number,...
    'include_long_profile',true,'title_flag',false);

d = 1;

fig_handle = figure(fig_number);

tile_handle = tiledlayout(fig_handle,2,1);
tile_handle.TileSpacing = 'none';
tile_handle.Padding = 'compact';
plot_position = [0.2857    0.0444    0.3979    0.3514];
fig_handle.Units = 'normalized';
fig_handle.OuterPosition = plot_position; %[0.1 0.3 0.8 0.5]
sigma_rexp = 0.536; %mm

hletter = 0.85;

%---------------- PLOT SIM 2D

timetemp = load([datadirs_sim{d},'densitytimeprofile.mat']);
simtimeprofile = timetemp.densitymatrix;
z_plot = linspace(timetemp.xlims(1),timetemp.xlims(2),size(simtimeprofile,2));
r_plot = linspace(0,timetemp.ylims(2),size(simtimeprofile,1))*10;


%first manual trimming
ind_r = (r_plot > trans_range(1)*10) & (r_plot < trans_range(2)*10);
density_plot_half = simtimeprofile(ind_r,:);
r_plot = r_plot(ind_r);

density_plot = [flip(density_plot_half);density_plot_half];
meanstd_density = 3*std(density_plot,[],'all');
max_opaqueness = 1;
ind_opaque = max_opaqueness*density_plot;
ind_opaque(density_plot > meanstd_density) = max_opaqueness*meanstd_density;
if d == 1
    max_indopaque_exp = max(ind_opaque,[],'all');
end
ind_opaque = ind_opaque/max_indopaque_exp;
if P.include_long_profile
    ax_density_sim(d) = nexttile;
else
    ax_density_sim(d) = axes('parent',fig_handle,'NextPlot','add','color','none');
end
% set limits to the axis
ax_density_sim(d).XLim = xi_range; %[min(P.SCxaxis),max(P.SCxaxis)];
ax_density_sim(d).YLim = [-max(r_plot),max(r_plot)];
ax_density_sim(d).XDir = 'reverse';
ax_density_sim(d).Box = 'on';

%     ax_density_sim(d).YTickMode = 'manual';


imagesc(ax_density_sim(d),'XData',z_plot,'YData',[-max(r_plot),max(r_plot)],'CData',double(density_plot>0),'alphadata',ind_opaque);
grad = colorGradient([1 1 1],[0 0 0],2);
colormap(ax_density_sim(d),grad);
text(0.025,hletter,letters{2*(d-1)+1},'Units','normalized','FontSize',12)
%         ylabel('x (mm)')

yline(sigma_rexp,'r','LineWidth',1)
yline(-sigma_rexp,'r','LineWidth',1)

ax_density_sim(d).XTickLabel = [];
%     ax_density_sim(d).XTick = [];
%         ax_density_sim(d).YTickLabel = [];
ax_density_sim(d).YTick = [-1,0,1];

ylabel('x (mm)');

%---------------- PLOT SIM LINEOUT

axlongprofile_sim(d) = nexttile;

% manual triming of experimental data
ind_r = (r_plot < sigma_rexp*10);
density_trim = density_plot(ind_r,:);
r_plot = r_plot(ind_r);

%end of manual triming

long_profile = sum(density_trim);
plongprofile = plot(z_plot,long_profile);
axlongprofile_sim(d).XDir = 'reverse';
xlim(xi_range); ylim([0 5e10])

ylabel('density (a.u.)');
xlabel(tile_handle,'delay (ps)');
drawnow;
pause(1);

P.fig_handle = fig_handle;
P.plot_name = 'g0density'; %'allprofiles';
P.save_plot();
%---------------------------------------------------------------------
% Plot waterfall, energy gain from Ez along xi (point-like bunch),
% maximum energy gain and position, and phase evolution.
%
% Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 21/09/2020
%---------------------------------------------------------------------
clear;

load('color_purple_to_green.mat');

%% input parameters
datadirs = {'gm20','gm15','gm10','gm5',...
    'g0','gp5','gp10','gp15','gp20'};

% datadirs = {'gp5','gp10'};

plasmaden = 1.81e14;
dump_list = 0:1:100;
useAvg = false;
dataformat = 'mat';
leg = {'-2 %/m','-1.5 %/m','-1 %/m','-0.5 %/m','0 %/m','0.5 %/m','1 %/m','1.5 %/m','2 %/m'};
xi_range = [5.7 6.3];
lineout_point = 5;
plots_dir = ['gradsim_paper/energy_gain'];
plot_name = ['energygain'];

OPA = OsirisPhaseAnalysis('datadir',datadirs{1},...
    'property','fields','wakefields_direction','long',...
    'plasmaden',plasmaden,...
    'dump_list',dump_list,'useAvg',useAvg,...
    'dataformat',dataformat,...
    'force_waterfall',false,...
    'lineout_point',lineout_point);
P = Plotty('plasmaden',plasmaden,'plot_name',plot_name,'plots_dir',plots_dir);


for d = 1:length(datadirs)
    OPA.datadir = datadirs{d};
    
    % waterfall plot
    OPA.build_waterfall();
    ind_xi = fliplr((OPA.waterfall_xi < xi_range(2)) & (OPA.waterfall_xi > xi_range(1)));
    waterfall_mat = OPA.waterfall_mat(:,ind_xi);
    waterfall_z = linspace(OPA.waterfall_z(1),OPA.waterfall_z(2),size(OPA.waterfall_mat,1));
    
    
    % field integral 
    field_integral = cumtrapz(waterfall_z,-(waterfall_mat));
      

    
    %% max energy gain 
    [max_energygain,pos_maxenergygain] = max(field_integral);
    
    [maxmax_energygain,best_position] = max(max_energygain);
    energy_gain(d,:) = field_integral(:,best_position);    
    field_amplitude(d,:) = -flip(waterfall_mat(:,best_position));

end

% energy gain
fig_eg = figure(1);
colororder(cc);
line_style = {'-','-','-','-','-','-.','-.','-.','-.'};

for d = 1:length(datadirs)
    hold on
    plot(waterfall_z,energy_gain(d,:),line_style{d},'LineWidth',2);
end
hold off

xlim([0 10])
ylim([-80 2000])

xlabel('distance from plasma end (m)')
ylabel ('energy gain (MV)')

legend(leg,'Location','best')
P.fig_handle = fig_eg;
P.save_plot();

% field amplitude
fig_fa = figure(2);
colororder(cc);
plot(waterfall_z,field_amplitude,'LineWidth',2);

xlim([0 10])

xlabel('z(m)')
ylabel ('E_z (MV)')

legend(leg,'Location','best')
P.fig_handle = fig_fa;
P.plot_name = 'field';
P.save_plot();


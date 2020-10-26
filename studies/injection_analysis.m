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
% Last update: 05/10/2020
%---------------------------------------------------------------------
clear;

% load('color_purple_to_green.mat');

%% input parameters
datadirs = {'R2gap0_2e14','R2gap100_2e14'};

plasmaden = 2e14;
dump_list = 0:1:200;
useAvg = true;
dataformat = 'h5';
leg = {'gap = 0 m','gap = 1 m'};
xi_range = [5.0 6.0];
dephasing_xi = 5.8;
lineout_point = 3;
plots_dir = ['gap/energy_gain'];
plot_name = ['energygain'];
xlim_plot = [0 20];

OPA = OsirisPhaseAnalysis('datadir',datadirs{1},...
    'property','fields','wakefields_direction','long',...
    'plasmaden',plasmaden,...
    'dump_list',dump_list,'useAvg',useAvg,...
    'dataformat',dataformat,...
    'force_waterfall',true,...
    'lineout_point',lineout_point,...
    'dephasing_xi',dephasing_xi,'dephasing_search','0x');
P = Plotty('plasmaden',plasmaden,'plot_name',plot_name,'plots_dir',plots_dir);

% initialization



for d = 1:length(datadirs)
    OPA.datadir = datadirs{d};
    switch OPA.datadir
        case 'R2gap0_2e14'
            OPA.dump_list = 0:1:200;
        case 'R2gap100_2e14'
            OPA.dump_list = 0:1:210;
    end
    
    % waterfall plot
    OPA.build_waterfall();
    ind_xi = fliplr((OPA.waterfall_xi < xi_range(2)) & (OPA.waterfall_xi > xi_range(1)));
    waterfall_mat = OPA.waterfall_mat(:,ind_xi);
    waterfall_z{d} = linspace(OPA.waterfall_z(1),OPA.waterfall_z(2),size(OPA.waterfall_mat,1));
    
    
    % field integral
    field_integral = cumtrapz(waterfall_z{d},-(waterfall_mat));
    
    
    
    % max energy gain
    [max_energygain,pos_maxenergygain] = max(field_integral);
    
    [maxmax_energygain,best_position] = max(max_energygain);
    energy_gain{d} = field_integral(:,best_position);
    
    % field at max energy gain
    field_amplitude{d} = -flip(waterfall_mat(:,best_position));
    
    % maximum fields
    [max_amplitude{d},ind_max_amplitude] = max(-flip(waterfall_mat),[],2);
    
    %% waterfall plot
    % waterfall plot datas
    P.waterfall_xi = OPA.waterfall_xi;
    P.waterfall_z = OPA.waterfall_z;
    P.waterfall_mat = OPA.waterfall_mat;
    P.property = OPA.property;
    
    P.fig_number = 10+d;
    P.waterfall_plot();
    x_range = xi_range;
    xlim(x_range);
    ind_x = fliplr((OPA.waterfall_xi < x_range(2)) & (OPA.waterfall_xi > x_range(1)));
    colormap_low = min(OPA.waterfall_mat(:,ind_x),[],'all');
    colormap_upp = max(OPA.waterfall_mat(:,ind_x),[],'all');
    caxis([colormap_low colormap_upp]);
    colormap(bluewhitered);
    ax = imgca(P.fig_handle);
   
    
        
    % phase
    dump_list_save =  OPA.dump_list;
    switch OPA.datadir
        case 'R2gap0_2e14'
            OPA.dump_list = [101:1:200];
        case 'R2gap100_2e14'
            OPA.dump_list = [111:1:210];
    end % switch datadir
    
    OPA.dephasing(); % out: dephasing_line
    dephasing_line(d,:) = OPA.dephasing_line;
    dephasing_z(d,:) = linspace(OPA.waterfall_z(1),OPA.waterfall_z(2),size(OPA.waterfall_mat,1));
    OPA.dump_list = dump_list_save;
    
    %% Add dephasing line to waterfall plot
    hold on
    plot((-OPA.dephasing_line)*OPA.plasma_wavelength+OPA.simulation_window-OPA.dephasing_first,...
        linspace(OPA.waterfall_z(1),OPA.waterfall_z(2),length(OPA.dephasing_line)),...
        'LineWidth',2,'Parent',ax,'color','k')
%     plot(bunches_x,bunches_y,'LineWidth',2,'Parent',ax,'color',[1, 0, 0])
    ax.XDir = 'reverse';
    ax.YDir = 'normal';
%     title(['\xi = ',num2str(dephasing_xi),' cm behind seed. pos. (',datadir,')'])
    hold off
    drawnow;
    
    P.plot_name = ['waterfalldephasingline',OPA.datadir];
    P.save_plot();
end




%% plot maximum energy gain
fig_eg = figure(1);
% colororder(cc);
for d = 1:length(datadirs)
    hold on
    plot(waterfall_z{d},energy_gain{d},'LineWidth',2);
    hold off
end

xlim(xlim_plot)
% ylim([-80 4000])

xlabel('distance from plasma end (m)')
ylabel ('energy gain (MeV)')

legend(leg,'Location','best')
P.fig_handle = fig_eg;
P.plot_name = 'energygain';
P.save_plot();


%% field at maximum energy gain
fig_fa = figure(2);
% colororder(cc);
for d = 1:length(datadirs)
    hold on
    plot(waterfall_z{d},field_amplitude{d},'LineWidth',2);
    hold off
end
xlim(xlim_plot)

xlabel('z (m)')
ylabel ('E_z (MV/m)')

legend(leg,'Location','best')
P.fig_handle = fig_fa;
P.plot_name = 'field';
P.save_plot();

%% maximum fields
fig_mf = figure(3);
% colororder(cc);
for d = 1:length(datadirs)
    hold on
    plot(waterfall_z{d},max_amplitude{d},'LineWidth',2);
    hold off
end
xlim(xlim_plot)

xlabel('z (m)')
ylabel ('max E_z (MV/m)')

legend(leg,'Location','best')
P.fig_handle = fig_mf;
P.plot_name = 'maxfield';
P.save_plot();

%% phase
fig_ph = figure(4);
% colororder(cc);
dephasing_z = linspace(0,10,length(dephasing_line));
plot(dephasing_z,dephasing_line','LineWidth',2);

xlim([0 10])
% ylim([-1 0])
xlabel('propagation (m)')
ylabel ('phase (\lambda_p)')

legend(leg,'Location','best')
P.fig_handle = fig_ph;
P.plot_name = 'dephasing';
P.save_plot();

%% phase 2
fig_ph = figure(5);
% colororder(cc);
dephasing_z = linspace(0,10,length(dephasing_line));
plot(dephasing_z,dephasing_line','LineWidth',2);

hold on
yline(0,'LineWidth',1)
yline(-0.25,'LineWidth',1)
hold off

xlim([0 10]);
ylim([-0.3 0.1])
xlabel('propagation (m)')
ylabel ('phase (\lambda_p)')

legend(leg,'Location','best','Autoupdate','off')
P.fig_handle = fig_ph;
P.plot_name = 'dephasing25';
P.save_plot();



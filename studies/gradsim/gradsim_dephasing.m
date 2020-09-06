%________________________________________________________________________
% Script to produce waterfall plots for different xi ranges of the proton bunch.
% Special version to produce the plot for the paper.
%
% Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 18/06/2020
%________________________________________________________________________
close all;

% datadirs = {'gm20','gm15','gm10','gm5','g0','gp5','gp10','gp15','gp20'};
datadirs = {'gm20','gm15','gm10','gm5','g0','gp1','gp2','gp5','gp10'};

% datadirs = {'gp2'};
plasmaden = 1.81e14;
dump_list = 0:100;
useAvg = true;
dataformat = 'h5';
dephasing_xi = 7; % cm

% load color order for 9 gradients
load('color_purple_to_green.mat');

for d = 1:length(datadirs)
    
    datadir = datadirs{d};
    
    switch datadir
        case  'gm20'
            title_g = 'g = -20 %';
        case 'gm15'
            title_g = 'g = -15 %';
        case 'gm10'
            title_g = 'g = -10 %';
        case 'gm5'
            title_g = 'g = -5 %';
        case  'g0'
            title_g = 'g = 0 %';
        case 'gp5'
            title_g = 'g = +5 %';
        case 'gp10'
            title_g = 'g = +10 %';
        case 'gp15'
            title_g = 'g = +15 %';
        case 'gp20'
            title_g = 'g = +20 %';
    end
    
    OPA = OsirisPhaseAnalysis('datadir',datadir,...
        'property','fields','wakefields_direction','long',...
        'plasmaden',plasmaden,...
        'dump_list',dump_list,'useAvg',useAvg,...
        'dataformat',dataformat,'dephasing_xi',dephasing_xi,...
        'force_waterfall',false);
    OPA.dephasing();
    fig1 = figure(1);
    colororder(cc)

    hold on
    p1 = plot(linspace(OPA.waterfall_z(1),OPA.waterfall_z(2),length(OPA.dephasing_line)),...
        OPA.dephasing_line,'LineWidth',2);
    hold off
    
    dephasing_z = linspace(OPA.waterfall_z(1),OPA.waterfall_z(2),...
        length(OPA.dephasing_line));
    
end
% P1 = Plotty();
% P1.plasmaden = linspace(1.81e14,(1-0.2)*1.81e14,length(OPA.dephasing_line));
% P1.den2freq();
% dephasing_frequency = P1.plasmafreq;
% 
% k = dephasing_frequency/P1.c_cm/(2*pi);
% phase_diff = k*7;
% 
% normalization = OPA.dephasing_line(50) - phase_diff(50)./OPA.plasma_wavelength;
% 
% phase_sin_plot = phase_diff/OPA.plasma_wavelength + normalization;
% hold on
% plot(dephasing_z,phase_sin_plot)
% hold off
% pause;


ylabel('dephasing (\lambda_p)')
xlabel('propagation distance (m)')
% legend('gm20','gm15','gm10','gm5','g0','gp5','gp10','gp15','gp20','location','best')
legend('gm20','gm15','gm10','gm5','g0','gp1','gp2','gp5','gp10','location','best')

xlim([0 10])
ylim([-7 1])
title(['\xi = ',num2str(dephasing_xi),' cm behind seed. pos.'])


P = Plotty('fig_handle',fig1);


P.plots_dir = 'dephasing/gradsim';
P.plot_name = ['dephasing',num2str(dephasing_xi),'x'];
% P.save_plot();
% 
% fig2 = figure(2);
% imagesc(OPA.waterfall_xi,OPA.waterfall_z,rot90(OPA.waterfall_mat,2));
% ax = imgca(fig2);
% ax.XDir = 'reverse';
% ax.YDir = 'normal';
% ylim(OPA.waterfall_z)
% xlim([6 8])
% colormap(bluewhitered)
% 
% hold on
% plot((-OPA.dephasing_line)*OPA.plasma_wavelength+OPA.simulation_window-OPA.dephasing_first,...
%     linspace(OPA.waterfall_z(1),OPA.waterfall_z(2),length(OPA.dephasing_line)),...
%     'LineWidth',2,'Parent',ax,'color','k')
% ax.XDir = 'reverse';
% ax.YDir = 'normal';
% hold off
% 
% xlabel('\xi (cm)')
% ylabel('propagation distance (m)')
% colorbar
% cbar.Label.String = 'E_z (MV/m)';
% P = Plotty('fig_handle',fig2);
% P.plots_dir = 'dephasing/gradsim';
% P.plot_name = ['dephasing',num2str(dephasing_xi)];
% P.save_plot();
% 
% 

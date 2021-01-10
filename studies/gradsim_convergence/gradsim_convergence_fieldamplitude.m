%________________________________________________________________________
% Calculate and plot the maximum fields
%
% Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 14/07/2020
%________________________________________________________________________

clear;
% close all;

datadirs = {'g0zh','g0rh','g0','g0z2','g0r2','g0dt92'};


load('color_purple_to_green.mat');
leg = {'-2 %/m','-1.5 %/m','-1 %/m','-0.5 %/m','0 %/m','0.5 %/m','1 %/m','1.5 %/m','2 %/m'};
% leg = {'0 %/m','0.5 %/m','1 %/m','1.5 %/m','2 %/m'};
% cc = cc(5:end,:);
% save file name
plots_dir = ['gradsim_convergence/amplitude/gradsim'];
saveplot = true;

% simulation parameters
plasmaden = 1.81e14; 
dump_list = 1:1:100;
sigma_z = 6.98516; % cm
bunch_center = 3.7;
wakefields_direction = 'trans';

search_type = 'mean';
search_xi = 6.98516;
trans_range = [0 100];%[0.002 0.01];
useAvg = false;
dataformat = 'mat';

% Initialization
amplitudes = cell(length(datadirs),1);
positions = cell(length(datadirs),1);
propagations = cell(length(datadirs),1);
for d = 1:length(datadirs)
    datadir = datadirs{d};
    switch datadirs{d}
        case 'g0'
            dataformat = 'mat';
        otherwise
            dataformat = 'h5';
    end
    OWA = OsirisWakefieldAnalysis('datadir',datadir,'dump_list',dump_list,...
       'wakefields_direction',wakefields_direction,...
        'useAvg',useAvg,...
        'search_type',search_type,...
        'plasmaden',plasmaden,'trans_range',trans_range,...
        'sigma_z',sigma_z,'bunch_center',bunch_center,...
        'dataformat',dataformat);
%     OWA.xi_range = search_xi + 0.6*[OWA.plasma_wavelength,-OWA.plasma_wavelength];
%     OWA.xi_range = [6.3 5.7];
    
    % find max field
    OWA.amplitude_vs_z();
    amplitudes{d} = OWA.denorm_Efield(OWA.amplitude_z);
    propagations{d} = OWA.propagation_z/100;
    positions{d} = OWA.pos_amplitude_z;
    
end

%% plot results 
P = Plotty('plots_dir',plots_dir);

fig_amplitude = figure(5);
line_style = {'-','-','-','-','-','-','-.','-.','-.'};
% colororder(cc);

hold on
for d = 1:length(datadirs)
    plot(propagations{d},(amplitudes{d}-amplitudes{3})./amplitudes{3}+1,line_style{d},'LineWidth',2);  
end
hold off
fig_amplitude.Units = 'normalized';
fig_amplitude.OuterPosition = [0 0.25 0.6 0.55];
axis('tight');
xlim([0 10]);
xlabel('z (m)');
% ylabel(['max. E_z (MV/m)']);
ylabel(['rel. mean defocusing fields (MV/m)']);
legend(datadirs,'location','best')
% P.plot_name = 'maxlongxi6';
P.plot_name = 'relmeandefocusing';
P.fig_handle = fig_amplitude;
P.save_plot();


fig_amplitude = figure;
colororder(cc)
hold on
for d = 1:length(datadirs)
    plot(propagations{d},positions{d},line_style{d},'LineWidth',2);  
end
hold off
fig2.Units = 'normalized';
fig2.OuterPosition = [0 0.25 0.6 0.55];
axis('tight');
xlim([0 10]);

% % title(['',meanormax,' E field position along bunch']);
% xlabel('z (cm)');
% ylabel(['max. E_z position ([z - z_c]/\sigma_z)']);
% legend(leg,'Location','best')
% P.plot_name = 'pos_max_long_xi7';
% P.fig_handle = fig2;
% P.save_plot();



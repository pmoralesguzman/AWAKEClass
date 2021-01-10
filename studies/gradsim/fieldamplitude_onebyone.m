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

% datadirs = {'gm15','gm10','gm5','g0','gp1','gp2','gp5','gp10','gp15'};
% datadirs = {'gm20','gm15','gm10','gm5','g0','gp5','gp10','gp15','gp20'};
datadirs = {'gp20'};
load('color_purple_to_green.mat');
% leg = {'-2 %/m','-1.5 %/m','-1 %/m','-0.5 %/m','0 %/m','0.5 %/m','1 %/m','1.5 %/m','2 %/m'};
leg = {'0 %/m','0.5 %/m','1 %/m','1.5 %/m','2 %/m'};
% cc = cc(5:end,:);
% save file name
plots_dir = ['gradsim/amplitude/'];
saveplot = true;
colornumber = 5;

% simulation parameters
plasmaden = 1.81e14;
dump_list = 0:1:100;
sigma_z = 6.98516;% 6.98516; % cm
bunch_center = 3.7;
wakefields_direction = 'long';

search_type = 'max';
search_xi = 6.98516;
trans_range = [0.003 0.005];
useAvg = false;
dataformat = 'mat';

% Initialization
amplitudes = cell(length(datadirs),1);
positions = cell(length(datadirs),1);
propagations = cell(length(datadirs),1);
for d = 1:length(datadirs)
    datadir = datadirs{d};
    
    OWA = OsirisWakefieldAnalysis('datadir',datadir,'dump_list',dump_list,...
        'wakefields_direction',wakefields_direction,...
        'useAvg',useAvg,...
        'search_type',search_type,...
        'plasmaden',plasmaden,'trans_range',trans_range,...
        'sigma_z',sigma_z,'bunch_center',bunch_center,...
        'dataformat',dataformat);
%         OWA.xi_range = search_xi + 0.6*[OWA.plasma_wavelength,-OWA.plasma_wavelength];
    OWA.xi_range = [18 0];
    
    % find max field
    OWA.amplitude_vs_z();
    amplitudes{d} = OWA.denorm_Efield(OWA.amplitude_z);
    propagations{d} = OWA.propagation_z/100;
    positions{d} = OWA.pos_amplitude_z;
    
end

%% plot results
P = Plotty('plots_dir',plots_dir);


for d = 1:length(datadirs)
    
    switch datadirs{d}
        case 'gp20'
            title_g = 'g = +2 %/m';
            grad_sim = 0.02;
        case 'gp15'
            title_g = 'g = +1.5 %/m';
            grad_sim = 0.015;
        case 'gp10'
            title_g = 'g = +1 %/m';
            grad_sim = 0.01;
        case 'gp5'
            title_g = 'g = +0.5 %/m';
            grad_sim = 0.005;
        case {'g0','g0d2'}
            title_g = 'g = 0 %/m';
            grad_sim = 0.0;
        case 'gm5'
            title_g = 'g = -0.5 %/m';
            grad_sim = -0.005;
        case 'gm10'
            title_g = 'g = -1 %/m';
            grad_sim = -0.01;
        case 'gm15'
            title_g = 'g = -1.5 %/m';
            grad_sim = -0.015;
        case 'gm20'
            title_g = 'g = -2 %/m';
            grad_sim = -0.02;
        case 'R2gap0_2e14'
            title_g = 'density step';
    end % switch datadir
    
    fig_amplitude = figure;
    if isempty(colornumber)
        colornumber = d;
    end
    plot(propagations{d},amplitudes{d},'LineWidth',2,'color',cc(colornumber,:));
    
    fig_amplitude.Units = 'normalized';
    fig_amplitude.OuterPosition = [0 0.25 0.4 0.4];
    axis('tight');
    xlim([0 10]);
    ylim([0 350]);
    xlabel('z (m)');
    ylabel(['max. E_z (MV/m) around \xi_0 = 7 cm']);
%     title(title_g);
    % ylabel(['mean defocusing fields (MV/m)']);
%     legend(leg,'Location','best')
    drawnow;
    P.plot_name = [datadirs{d},'maxlongxi7_g0'];
    P.fig_handle = fig_amplitude;
    P.save_plot();
end

for d = 1:length(datadirs)
    fig_position = figure;
    plot(propagations{d},positions{d},'LineWidth',2);

fig_position.Units = 'normalized';
fig_position.OuterPosition = [0 0.25 0.6 0.55];
axis('tight');
xlim([0 10]);

% % title(['',meanormax,' E field position along bunch']);
xlabel('z (cm)');
ylabel(['max. E_z position ([z - z_c]/\sigma_z)']);
legend(leg,'Location','best')
drawnow;

P.plot_name = 'pos_max_long_xi7';
P.fig_handle = fig_position;
P.save_plot();
end


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
close all;

% datadirs = {'gm20','gm15','gm10','gm5','g0','gp5','gp10','gp15','gp20'};
datadirs = {'gm15','gm10','gm5','g0','gp1','gp2','gp5','gp10','gp15'};
grads_sim = [-20,-15,-10,-5,0,5,10,15,20];

% save file name
plots_dir = ['amplitude/gradsim'];
saveplot = false;

% simulation parameters
plasmaden = 1.81e14; 
dump_list = 1:1:100;
sigma_z = 6.98516; % cm
bunch_center = 3.7;
wakefields_direction = 'long';

search_type = 'max';
search_xi = 6.98516;
trans_range = [0.004 0.2];
useAvg = true;
dataformat = 'h5';

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
%     OWA.xi_range = search_xi + [OWA.plasma_wavelength/2,-OWA.plasma_wavelength/2];
    OWA.xi_range = [21 0];
    
    % find max field
    OWA.amplitude_vs_z();
    amplitudes{d} = OWA.denorm_Efield(OWA.amplitude_z);
    propagations{d} = OWA.propagation_z/100;
    positions{d} = OWA.pos_amplitude_z;
    
end

%% plot results 
P = Plotty('plots_dir',plots_dir);

% color_order = get(gca,'colororder');

cc = [118,42,131
153,112,171
194,165,207
231,212,232
10,10,10
217,240,211
166,219,160
90,174,97
27,120,55]/256;

fig = figure(1);
colororder(cc);
% set(0,'DefaultAxesColorOrder',brewermap(9,'Blues')) 


hold on
for d = 1:length(datadirs)
    plot(propagations{d},amplitudes{d},'LineWidth',2);  
end
hold off
fig.Units = 'normalized';
fig.OuterPosition = [0 0.25 0.6 0.55];
axis('tight');
xlim([0 10]);
xlabel('propagation distance (m)');
% ylabel(['max. E_z (E_0) in the first 2 \sigma']);
ylabel(['mean defocusing fields (MV/m)']);
% legend('gm20','gm15','gm10','gm5','g0','gp5','gp10','gp15','gp20','Location','best')
legend('gm15','gm10','gm5','g0','gp1','gp2','gp5','gp10','gp15','Location','best')
P.plot_name = 'mean_trans_square';
P.fig_handle = fig;
P.save_plot();
% ylim([0 250]);
% annotation('textarrow',[.25 .28],[.7 .65],'LineWidth',1.5,'HeadStyle','none','String','3.7 nC','Color',color_order(5,1:3))
% annotation('textarrow',[.3 .33],[.7 .65],'LineWidth',1.5,'HeadStyle','none','String','921 pC','Color',color_order(4,1:3))
% annotation('textarrow',[.35 .38],[.7 .65],'LineWidth',1.5,'HeadStyle','none','String','184 pC','Color',color_order(3,1:3))
% annotation('textarrow',[.41 .44],[.7 .65],'LineWidth',1.5,'HeadStyle','none','String','cut bunch','Color',color_order(1,1:3))
% annotation('textarrow',[.47 .5],[.7 .65],'LineWidth',1.5,'HeadStyle','none','String','no seed','Color',color_order(2,1:3))


% title(['R2 no density step (full window)']);



% leg = {'ionization front (13.5 MV/m)', 'no seed bunch (0.15 MV/m)',...
%     'seed bunch charge \times 1 (12.5 MV/m)', 'seed bunch \times 5 (65.2 MV/m)',...
%     'seed bunch charge \times 20 (298.5 MV/m)'};
% seeding long. fields: 
% ionization front: 13.5 MV/m, 
% no seed bunch: 0.15 MV/m
% seed bunch x1: 12.5 MV/m
% seed bunch x5: 65.2 MV/m 
% seed bunch x20: 298.5 MV/m 


% if ~isfolder(plots_dir); mkdir(plots_dir); end
% if saveplot
%     saveas(gcf,[plots_dir,'maxfieldR2nods_full.png'])
%     savefig([plots_dir,'maxfieldR2nods_full.fig'])
% end

fig2 = figure(2);
colororder(cc)
hold on
for d = 1:length(datadirs)
    plot(propagations{d},positions{d},'LineWidth',2);  
end
hold off
fig2.Units = 'normalized';
fig2.OuterPosition = [0 0.25 0.6 0.55];
axis('tight');
xlim([0 10]);

% % title(['',meanormax,' E field position along bunch']);
xlabel('propagation distance (cm)');
ylabel(['max. E_z position ([z - z_c]/\sigma_z)']);
% legend('gm20','gm15','gm10','gm5','g0','gp5','gp10','gp15','gp20')
legend('gm15','gm10','gm5','g0','gp1','gp2','gp5','gp10','gp15')
P.plot_name = 'pos_max_xi7';
P.fig_handle = fig2;
P.save_plot();
% 
% % leg = {'ionization front (13.5 MV/m)', 'no seed bunch (0.15 MV/m)',...
% %     'seed bunch charge \times 1 (12.5 MV/m)', 'seed bunch \times 5 (65.2 MV/m)',...
% %     'seed bunch charge \times 20 (298.5 MV/m)'};legend(leg,'location','best');
% 
% if saveplot
%     saveas(gcf,[plots_dir,'maxfield_pos_R2nods.png'])
%     savefig([plots_dir,'maxfield_pos_R2nods.fig'])
% end

% figure
% bar(grads_sim,total_field);


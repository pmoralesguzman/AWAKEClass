%________________________________________________________________________
% Calculate and plot the mean transverse fields, together with manually
% positioned tags.
%
% Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 11/01/2021
%________________________________________________________________________

answer = questdlg('Are you sure?');
switch answer
    case 'Yes'
        % nothing
    otherwise
        return;
end

clear;
% close all;

run_it = 1;

% datadirs = {'gm15','gm10','gm5','g0','gp1','gp1,'gp5','gp10','gp15'};
datadirs = {'gm20','gm15','gm10','gm5','g0','gp5','gp10','gp15','gp20'};
% datadirs = {'g0'};

load('color_purple_to_green.mat');
load('color_red_to_blue.mat');

leg = {'-2 %/m','-1.5 %/m','-1 %/m','-0.5 %/m','0 %/m','0.5 %/m','1 %/m','1.5 %/m','2 %/m'};
% leg = {'0 %/m','0.5 %/m','1 %/m','1.5 %/m','2 %/m'};
% ccrb = ccrb(5:end,:);
% save file name
plots_dir = ['gradsim_paper/amplitude/all'];
saveplot = true;

% simulation parameters
plasmaden = 1.81e14;
dump_list = 0:1:100;
wakefields_direction = 'trans';

search_type = 'mean';
trans_range = [0 100];%[0.002 0.01];
useAvg = false;
dataformat = 'mat';

% Initialization
amplitudes = cell(length(datadirs),1);
positions = cell(length(datadirs),1);
propagations = cell(length(datadirs),1);

if run_it ==1
    
    for d = 1:length(datadirs)
        datadir = datadirs{d};
        
        OWA = OsirisWakefieldAnalysis('datadir',datadir,'dump_list',dump_list,...
            'wakefields_direction',wakefields_direction,...
            'useAvg',useAvg,...
            'search_type',search_type,...
            'plasmaden',plasmaden,'trans_range',trans_range,...
            'dataformat',dataformat);
        
        
        % find max field
        OWA.amplitude_vs_z();
        amplitudes{d} = OWA.denorm_Efield(OWA.amplitude_z);
        propagations{d} = OWA.propagation_z/100;
        positions{d} = OWA.pos_amplitude_z;
        
    end
    
else
    load('loading_files/gradsim_fieldamplitude_transmean.mat');
end % if run it

%% plot results
P = Plotty('plots_dir',plots_dir);
fontsize_annotation = 12;
fontsize_label = 14;

fig_amplitude = figure(5);
colororder(ccrb);
line_style = {':','--','-.','-','-','-','-.','--',':'};

hold on
for d = 1:length(datadirs)
    plot(propagations{d},amplitudes{d},line_style{d},'LineWidth',2);
end
hold off
ax = gca;
ax.FontSize = fontsize_label;
fig_amplitude.Units = 'normalized';
fig_amplitude.OuterPosition = [0 0.25 0.6 0.55];
axis('tight');
xlim([0 10]);
xlabel('z (m)','FontSize',fontsize_label);
% ylabel(['max. E_z (MV/m)']);
ylabel(['mean defocusing fields (MV/m)'],'FontSize',fontsize_label);
addtextlines = true;
if addtextlines
    % Create textarrow
    annotation(fig_amplitude,'textarrow',[0.514499121265378 0.484182776801406],...
        [0.857427716849452 0.835493519441675],'String','g = 0 \%/m',...
        'Interpreter','latex',...
        'LineWidth',1,'Color',ccrb(5,:),...
        'HeadStyle','none','FontSize',fontsize_annotation);
    
    % Create textarrow
    annotation(fig_amplitude,'textarrow',[0.847539543057997 0.864235500878735],...
        [0.712861415752742 0.76271186440678],'String','g = 0.5 \%/m',...
        'Interpreter','latex',...
        'LineWidth',1,'Color',ccrb(6,:),...
        'HeadStyle','none','FontSize',fontsize_annotation);
    
    % Create textarrow
    annotation(fig_amplitude,'textarrow',[0.651362038664323 0.621045694200351],...
        [0.755234297108675 0.733300099700898],'String','g = 1 \%/m',...
        'Interpreter','latex',...
        'LineWidth',1,'Color',ccrb(7,:),...
        'HeadStyle','none','FontSize',fontsize_annotation);
    
    % Create textarrow
    annotation(fig_amplitude,'textarrow',[0.544346221441125,0.553163444639719],...
        [0.633190428713858,0.542372881355932],'String','g = 1.5 \%/m',...
        'Interpreter','latex',...
        'LineWidth',1,'Color',ccrb(8,:),...
        'HeadStyle','none','FontSize',fontsize_annotation);
    
    % Create textarrow
    annotation(fig_amplitude,'textarrow',[0.490553602811951 0.510544815465729],...
        [0.32751744765703 0.383848454636092],'String','g = 2 \%/m',...
        'Interpreter','latex',...
        'LineWidth',1,'Color',ccrb(9,:),...
        'HeadStyle','none','FontSize',fontsize_annotation);
    
    % Create textarrow
    annotation(fig_amplitude,'textarrow',[0.320738137082601 0.371265377855887],...
        [0.780658025922233 0.684945164506481],'String','g = -1 \%/m',...
        'Interpreter','latex',...
        'LineWidth',1,'Color',ccrb(3,:),...
        'HeadStyle','none','FontSize',fontsize_annotation);
    
    % Create textarrow
    annotation(fig_amplitude,'textarrow',[0.350175746924429 0.397627416520211],...
        [0.860418743768694 0.751744765702891],'String','g = -0.5 \%/m',...
        'Interpreter','latex',...
        'LineWidth',1,'Color',ccrb(4,:),...
        'HeadStyle','none','FontSize',fontsize_annotation);
    
    % Create textarrow
    annotation(fig_amplitude,'textarrow',[0.300087873462214 0.355448154657293],...
        [0.708873379860419 0.612163509471585],'String','g = -1.5 \%/m',...
        'Interpreter','latex',...
        'LineWidth',1,'Color',ccrb(2,:),...
        'HeadStyle','none','FontSize',fontsize_annotation);
    
    % Create textarrow
    annotation(fig_amplitude,'textarrow',[0.287785588752197 0.349736379613357],...
        [0.634097706879362 0.542372881355932],'String','g = -2 \%/m',...
        'Interpreter','latex',...
        'LineWidth',1,'Color',ccrb(1,:),...
        'HeadStyle','none','FontSize',fontsize_annotation);
end % if add text lines

% legend(leg,'Location','best')
P.plot_name = 'meandefocusing';
P.fig_handle = fig_amplitude;
P.save_plot();

save('loading_files/gradsim_fieldamplitude_transmean.mat','propagations','amplitudes');




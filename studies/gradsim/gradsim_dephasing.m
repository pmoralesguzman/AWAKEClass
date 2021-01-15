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
clear;
close all;

datadirs = {'gm20','gm15','gm10','gm5','g0','gp5','gp10','gp15','gp20'};
% datadirs = {'gm20','gm15','gm10','gm5','g0','gp1','gp2','gp5','gp10'};
leg = {'-2 %/m','-1.5 %/m','-1 %/m','-0.5 %/m','0 %/m','0.5 %/m','1 %/m','1.5 %/m','2 %/m'};

% datadirs = {'gp2'};
plasmaden = 1.81e14;
dump_list = 0:1:100;
useAvg = false;
dataformat = 'mat';
dephasing_xi = [15,7,1]; % cm

plots_dir = ['gradsim_paper/dephasing_lines/'];
plot_name = ['dephasing','1471','allgrads'];
P = Plotty('plots_dir',plots_dir,'plot_name',plot_name);


% load color order for 9 gradients
load('color_purple_to_green.mat');

OPA = OsirisPhaseAnalysis('datadir',datadirs{1},...
    'property','fields','wakefields_direction','long',...
    'plasmaden',plasmaden,...
    'dump_list',dump_list,'useAvg',useAvg,...
    'dataformat',dataformat,'dephasing_xi',dephasing_xi(1),...
    'force_waterfall',false); OPA.dephasing();

dephasing_z = zeros(length(datadirs),length(OPA.dephasing_line));
dephasing_lines = dephasing_z;
%% 
fig_dephase = figure(1);
fig_dephase.OuterPosition = [100 100 1200 400];
colororder(cc);
tt = tiledlayout(1,3);
tt.TileSpacing = 'compact';
tt.Padding = 'compact';

for xi = 1:length(dephasing_xi)
    
    for d = 1:length(datadirs)
        
        OPA = OsirisPhaseAnalysis('datadir',datadirs{d},...
            'property','fields','wakefields_direction','long',...
            'plasmaden',plasmaden,...
            'dump_list',dump_list,'useAvg',useAvg,...
            'dataformat',dataformat,'dephasing_xi',dephasing_xi(xi),...
            'force_waterfall',false);
        
        OPA.dephasing();
        
        
        dephasing_z(d,:) = linspace(OPA.waterfall_z(1),OPA.waterfall_z(2),...
            length(OPA.dephasing_line));
        dephasing_lines(d,:) = OPA.dephasing_line;
        OPA.progress_dump('dephasing lines',d,length(datadirs))
        
    end
    
    ax(xi) = nexttile;
    line_style = {'-','-','-','-','-','-.','-.','-.','-.'};
    hold on
    for d = 1:length(datadirs)
        p1 = plot(dephasing_z(d,:),dephasing_lines(d,:),line_style{d},'LineWidth',2);
    end
    yline(-1,'--','LineWidth',1,'color',[0 0.4470 0.7410])
    xline(4,'--','LineWidth',1,'color',[0 0.4470 0.7410]);
    hold off
    
    switch dephasing_xi(xi)
        case 15
            position_word = '(back of the bunch)';
        case 7
            position_word = '(middle of the bunch)';
        case 1
            position_word = '(front of the bunch)';
    end
    title(['\xi_0 = ',num2str(dephasing_xi(xi)),' cm ',position_word]);
    
    xlim([0 10])
    ylim([-6.5 0.5])
    
end

legend(ax(3),leg,'location','southeast','FontSize',10)
ylabel(tt,'zero-crossing position shift (\lambda_p)')
xlabel(tt,'z (m)')


P.fig_handle = fig_dephase;
P.save_plot();


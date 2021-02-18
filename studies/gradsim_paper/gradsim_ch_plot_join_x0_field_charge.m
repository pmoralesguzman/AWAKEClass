%________________________________________________________________________
% SCript joining the dephasing, charge and mean transverse field amplitude
% for a 3x3 tile plot with all info.
%
% Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 02/02/2020
%________________________________________________________________________

clear;
close all;

datadirs = {'gm20','gm15','gm10','gm5','g0','gp5','gp10','gp15','gp20'};
leg = {'-2 %/m','-1.5 %/m','-1 %/m','-0.5 %/m','0 %/m','0.5 %/m','1 %/m','1.5 %/m','2 %/m'};

plasmaden = 1.81e14;
dump_list = 0:1:100;
useAvg = false;
dataformat = 'mat';
dephasing_xi = [14,7,1]; % cm

plots_dir = ['gradsim_test/join/'];
plot_name = ['chargefielddephasing','xi1471','allgrads0'];
P = Plotty('plots_dir',plots_dir,'plot_name',plot_name);

% load color order for 9 gradients
load('color_red_to_blue.mat');
hletter = 0.85;
i_letter = 1;

letters = {'a)','b)','c)',...
    'd)','e)','f)','g)','h)','i)','j)','k)','l)'};

load('loading_files/gradsim_dephasing.mat')


OPA = OsirisPhaseAnalysis('datadir',datadirs{1},...
    'property','fields','wakefields_direction','long',...
    'plasmaden',plasmaden,...
    'dump_list',dump_list,'useAvg',useAvg,...
    'dataformat',dataformat,'dephasing_xi',dephasing_xi(1),...
    'force_waterfall',false); OPA.dephasing();

%%
fig_dephase = figure(1);
fig_dephase.OuterPosition = [100 100 1200 800];
colororder(ccrb);
tt = tiledlayout(3,3);
tt.TileSpacing = 'compact';
tt.Padding = 'compact';

fontsize_annotation = 12;
fontsize_label = 14;

for xi = 1:length(dephasing_xi)
    
    axi(xi) = nexttile;
    axi(xi).FontSize = fontsize_label;
    
    line_style = {':','--','-.','-','-','-','-.','--',':'};
    hold on
    for d = 1:length(datadirs)
        p1 = plot(squeeze(dephasing_z(xi,d,:)),squeeze(dephasing_lines(xi,d,:)),line_style{d},'LineWidth',2);
    end
    text(0.03,0.13,letters{i_letter},'Units','normalized','FontSize',fontsize_label)
    i_letter = i_letter + 1;
    yline(-1,'--','LineWidth',1,'color',[0 0 0])
    yline(-2,'--','LineWidth',1,'color',[0 0 0])
    xline(4,'--','LineWidth',1,'color',[0 0 0]);
    hold off
    
    switch dephasing_xi(xi)
        case 14
            position_word = '(back of the bunch)';
        case 7
            position_word = '(middle of the bunch)';
        case 1
            position_word = '(front of the bunch)';
    end % dephasing xi (title)
    title(['\xi_0 = ',num2str(dephasing_xi(xi)),' cm ',position_word]);
    
    xlim([0 10])
    
    if OPA.dephasing_first == 40
        ylim([-6.5 1.5]*2)
    else
        ylim([-4 0.5]*2)
    end % if OPA dephasing first
    
end % xi

legend(axi(3),leg,'location','southeast','FontSize',fontsize_annotation,'NumColumns',2)
ylabel(axi(1),{'zero-crossing','position shift (\lambda_p/2)'},'FontSize',fontsize_label)

axi(2).YTickLabel = [];
axi(3).YTickLabel = [];
axi(1).XTickLabel = [];
axi(2).XTickLabel = [];
axi(3).XTickLabel = [];


if OPA.dephasing_first == 40
    
else
    
% Create textarrow
annotation(fig_dephase,'textarrow',[0.530314371257486 0.593562874251498],...
    [0.727554179566565 0.773993808049537],...
    'Color',[0.79296875 0.09375 0.11328125],...
    'String',{'negative','gradient values'},...
    'Interpreter','latex',...
    'HeadStyle','none',...
    'FontSize',12);

% Create textarrow
annotation(fig_dephase,'textarrow',[0.471182634730539 0.566616766467066],...
    [0.759133126934985 0.852012383900929],'String','g = 0 \%/m',...
    'Interpreter','latex',...
    'HeadStyle','none',...
    'FontSize',12);

% Create textarrow
annotation(fig_dephase,'textarrow',[0.478667664670659 0.546781437125749],...
    [0.820433436532508 0.888544891640871],...
    'Color',[0.12890625 0.44140625 0.70703125],...
    'String',{'positive','gradient values'},...
    'Interpreter','latex',...
    'HeadStyle','none',...
    'FontSize',12);
end

%----------------------------------------------------------------------
%----------------------------------------------------------------------
%----------------------------------------------------------------------
%----------------------------------------------------------------------
%----------------------------------------------------------------------
load('loading_files/gradsim_fieldamplitude_transmean.mat')

%% plotting mean trans field

for x = 1:length(dephasing_xi)
    axa(x) = nexttile;
    axa(x).FontSize = fontsize_label;
    for d = 1:length(datadirs)
        hold on
        plot(squeeze(plot_z(x,d,:)),squeeze(fieldvsz(x,d,:)),...
            line_style{d},'LineWidth',2,'color',ccrb(d,:))
        hold off
    end % datadir
    text(0.03,0.9,letters{i_letter},'Units','normalized','FontSize',fontsize_label)
    i_letter = i_letter + 1;
    hold on
    xline(4,'--','LineWidth',1,'color',[0 0 0]);
    hold off
    
    ylim([0,max(fieldvsz,[],'all')])
    xlim([0,10])
    
end % xi

ylabel(axa(1),{'mean defocusing','field (MV/m)'})



axa(2).YTickLabel = [];
axa(3).YTickLabel = [];
axa(1).XTickLabel = [];
axa(2).XTickLabel = [];
axa(3).XTickLabel = [];

%----------------------------------------------------------------------
%----------------------------------------------------------------------
%----------------------------------------------------------------------
%----------------------------------------------------------------------
%----------------------------------------------------------------------

load('loading_files/gradsim_chargeinxi.mat');

i_color = [1:9];

%% plotting charge

for x = 1:length(dephasing_xi)
    ax(x) = nexttile;
    ax(x).FontSize = fontsize_label;
    charge_norm_xi = mean(squeeze(chargevsz(x,:,1)));
    for d = 1:length(datadirs)
        hold on
        plot(squeeze(plot_z(x,d,:)),squeeze(chargevsz(x,d,:)./charge_norm_xi),...
            line_style{d},'LineWidth',2,'color',ccrb(i_color(d),:))
        hold off
    end % datadir
    text(0.03,0.13,letters{i_letter},'Units','normalized','FontSize',fontsize_label)
    i_letter = i_letter + 1;
    hold on
    xline(4,'--','LineWidth',1,'color',[0 0 0]);
    hold off
    
    ylim([0,1.02])
    xlim([0,10])
    
end % xi

ylabel(ax(1),'charge fraction');

xlabel(tt,'z (m)')

ax(2).YTickLabel = [];
ax(3).YTickLabel = [];

P.fig_handle = fig_dephase;
P.save_plot();










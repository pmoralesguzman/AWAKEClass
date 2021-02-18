%________________________________________________________________________
% gradsim paper
% Script joining the dephasing, charge and mean transverse field amplitude
% for a 3x3 tile plot with all info. It does not calculate anything, bith
% rather takes the data from
% - gradsim_ch_mean_trans_fields_in_xi
% - gradsim_ch_charge_fraction_in_xi
% - gradsim_ch_x0shift_in_xi
%
% Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 18/02/2020
%________________________________________________________________________

clear;
close all;

plots_dir = ['gradsim_paper/join/'];
plot_name = ['chargefielddephasing','xi1471','allgrads0'];

% load color order for 9 gradients
load('color_red_to_blue.mat'); % loaded vars: ccrb

% load x0 shift data
load('loading_files/gradsim_dephasing.mat'); % loaded vars: dephasing_z, dephasing_lines

% plotting parameters
fontsize_annotation = 12;
fontsize_label = 14;
dephasing_first = 0;
letterbox_x = 0.03;
letterbox_y = 0.13;
letterbox_y2 = 0.09;

% cell plotting parameters
datadirs = {'gm20','gm15','gm10','gm5','g0','gp5','gp10','gp15','gp20'};
leg = {'-2 %/m','-1.5 %/m','-1 %/m','-0.5 %/m','0 %/m','0.5 %/m','1 %/m','1.5 %/m','2 %/m'};
letters = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)'};
line_style = {':','--','-.','-','-','-','-.','--',':'};

% study parameters
dephasing_xi = [14,7,1]; % cm

% initialize counters
i_letter = 1;

% initialize classes 
P = Plotty('plasmaden',1.81e14,'plots_dir',plots_dir,'plot_name',plot_name);

%% start the plotting, with dephasing
fig_join = figure(1);
fig_join.OuterPosition = [100 100 1200 800];
colororder(ccrb);
tt = tiledlayout(3,3);
tt.TileSpacing = 'compact';
tt.Padding = 'compact';

% dephasing (x0 shift)

for xi = 1:length(dephasing_xi)
    
    ax_x0(xi) = nexttile; 
    ax_x0(xi).FontSize = fontsize_label;
    
    hold on
    for d = 1:length(datadirs)
        p1 = plot(squeeze(dephasing_z(xi,d,:)),squeeze(dephasing_lines(xi,d,:)),...
            line_style{d},'LineWidth',2);
    end
    text(letterbox_x,letterbox_y,letters{i_letter},'FontSize',fontsize_label,'Units','normalized')
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
    ylim([-4 0.5]*2)
    
end % xi

legend(ax_x0(3),leg,'location','southeast','FontSize',fontsize_annotation,'NumColumns',2)
ylabel(ax_x0(1),{'zero-crossing','position shift (\lambda_p/2)'},'FontSize',fontsize_label)

ax_x0(2).YTickLabel = [];
ax_x0(3).YTickLabel = [];
ax_x0(1).XTickLabel = [];
ax_x0(2).XTickLabel = [];
ax_x0(3).XTickLabel = [];

% Create textarrow
annotation(fig_join,'textarrow',[0.530314371257486 0.593562874251498],...
    [0.727554179566565 0.773993808049537],...
    'Color',[0.79296875 0.09375 0.11328125],...
    'String',{'negative','gradient values'},...
    'Interpreter','latex',...
    'HeadStyle','none',...
    'FontSize',12);

% Create textarrow
annotation(fig_join,'textarrow',[0.471182634730539 0.566616766467066],...
    [0.759133126934985 0.852012383900929],'String','g = 0 \%/m',...
    'Interpreter','latex',...
    'HeadStyle','none',...
    'FontSize',12);

% Create textarrow
annotation(fig_join,'textarrow',[0.478667664670659 0.546781437125749],...
    [0.820433436532508 0.888544891640871],...
    'Color',[0.12890625 0.44140625 0.70703125],...
    'String',{'positive','gradient values'},...
    'Interpreter','latex',...
    'HeadStyle','none',...
    'FontSize',12);

%----------------------------------------------------------------------
%----------------------------------------------------------------------
%----------------------------------------------------------------------

%% mean trans field
load('loading_files/gradsim_fieldamplitude_transmean.mat'); % loaded var: fieldvsz, plot_z

for xi = 1:length(dephasing_xi)
    ax_fld(xi) = nexttile;
    ax_fld(xi).FontSize = fontsize_label;
    
    hold on
    for d = 1:length(datadirs)
        plot(squeeze(plot_z(xi,d,:)),squeeze(fieldvsz(xi,d,:)),...
            line_style{d},'LineWidth',2) % 'color',ccrb(d,:)
    end % datadir
    hold off
    
    text(letterbox_x,letterbox_y2,letters{i_letter},'FontSize',fontsize_label,'Units','normalized')
    i_letter = i_letter + 1;
    
    hold on
    xline(4,'--','LineWidth',1,'color',[0 0 0]);
    hold off
    
    ylim([0,max(fieldvsz,[],'all')])
    xlim([0,10])
    
end % xi

ylabel(ax_fld(1),{'mean defocusing','field (MV/m)'})

ax_fld(2).YTickLabel = [];
ax_fld(3).YTickLabel = [];
ax_fld(1).XTickLabel = [];
ax_fld(2).XTickLabel = [];
ax_fld(3).XTickLabel = [];

%----------------------------------------------------------------------
%----------------------------------------------------------------------
%----------------------------------------------------------------------

%% plotting charge
load('loading_files/gradsim_chargeinxi.mat'); % loaded vars: plot_z, chargevsz

for xi = 1:length(dephasing_xi)
    ax_ch(xi) = nexttile;
    ax_ch(xi).FontSize = fontsize_label;
    charge_norm_xi = mean(squeeze(chargevsz(xi,:,1)));
    
    for d = 1:length(datadirs)
        hold on
        plot(squeeze(plot_z(xi,d,:)),squeeze(chargevsz(xi,d,:)./charge_norm_xi),...
            line_style{d},'LineWidth',2); % 'color',ccrb(i_color(d),:)
        hold off
    end % datadir
    text(letterbox_x,letterbox_y,letters{i_letter},'FontSize',fontsize_label,'Units','normalized')
    i_letter = i_letter + 1;
    hold on
    xline(4,'--','LineWidth',1,'color',[0 0 0]);
    hold off
    
    ylim([0,1.02])
    xlim([0,10])
    
end % xi

ylabel(ax_ch(1),'charge fraction');
xlabel(tt,'z (m)')

ax_ch(2).YTickLabel = [];
ax_ch(3).YTickLabel = [];

P.fig_handle = fig_join;
P.save_plot();










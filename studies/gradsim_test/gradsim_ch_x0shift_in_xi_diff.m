%________________________________________________________________________
% gradsim
% Script to follow the zero crossing of fields from waterfall plots. It
% also creates the waterfall plots if not yet made. It also calculates the
% difference between xi0 lines at each z.
%
% Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 19/02/2021
%________________________________________________________________________
clear;
close all;

plots_dir = ['gradsim_test/x0shift_diff/'];

% load color order for 9 gradients
load('color_purple_to_green.mat'); % cc
load('color_red_to_blue.mat'); % ccrb

% cell plotting parameters
datadirs = {'gm20','gm15','gm10','gm5','g0','gp5','gp10','gp15','gp20'};
% datadirs = {'gm20','gm15','gm10','gm5','g0','gp1','gp2','gp5','gp10'};
% datadirs = {'gp5'};
leg = {'-2 %/m','-1.5 %/m','-1 %/m','-0.5 %/m','0 %/m','0.5 %/m','1 %/m','1.5 %/m','2 %/m'};
line_style = {':','--','-.','-','-','-','-.','--',':'};
colors = {'r','k',[0 0.5 0],'b',[0.5,0,0.5]};

% plotting parameters
fontsize_annotation = 12;
fontsize_label = 14;

% study parameters
plasmaden = 1.81e14;
dump_list = 0:1:100;
prop_dis = 3:1:100*9.958450708031729;
useAvg = false;
dataformat = 'mat';
dephasing_xi_ini = [1]; % cm

% initialize classes
P = Plotty('plots_dir',plots_dir);
OPA = OsirisPhaseAnalysis('datadir',datadirs{1},...
    'property','fields','wakefields_direction','long',...
    'plasmaden',plasmaden,...
    'dump_list',dump_list,'useAvg',useAvg,...
    'dataformat',dataformat,'dephasing_xi',dephasing_xi_ini,...
    'force_waterfall',false);

dephasing_xi = [0.9:OPA.plasma_wavelength:15]; % cm
dephasing_xi0_plot = dephasing_xi(1:end-1) + OPA.plasma_wavelength/2;
dephasing_lines_number = 1:(length(dephasing_xi) - 1);

OPA.dephasing();

% initialize variables
dephasing_z = zeros(length(OPA.dephasing_line),1);
dephasing_lines = zeros(length(dephasing_xi),length(datadirs),length(OPA.dephasing_line));
% dephasing_lines = zeros(length(dephasing_xi),length(OPA.dephasing_line));

%% start the script
% fig_dephase = figure(1);
% fig_dephase.OuterPosition = [100 100 1200 400];
% colororder(ccrb);
% tt = tiledlayout(1,3);
% tt.TileSpacing = 'compact';
% tt.Padding = 'compact';

for d = 1:length(datadirs)
    
    OPA.datadir = datadirs{d};
    for xi = 1:length(dephasing_xi)
        OPA.dephasing_xi = dephasing_xi(xi);
        OPA.dephasing_first = [];
        
        OPA.dephasing();
        if sum(dephasing_z) == 0
            dephasing_z = linspace(OPA.waterfall_z(1),OPA.waterfall_z(2),...
                length(OPA.dephasing_line));
        end % sum dephasing z
        
        dephasing_lines(xi,d,:) = -OPA.dephasing_line*2 + (OPA.simulation_window - OPA.dephasing_first)*2/OPA.plasma_wavelength;% HALF PLASMA WAVELENGTH !!!
%         dephasing_lines(xi,:) = -OPA.dephasing_line*2 + (OPA.simulation_window - OPA.dephasing_first)*2/OPA.plasma_wavelength; % HALF PLASMA WAVELENGTH !!!

        OPA.progress_dump('dephasing lines',xi,length(dephasing_xi))

    end % dephasing
     OPA.progress_dump('DATADIR',d,length(datadirs))
end % datadirs

dephasing_diffs_all = diff(dephasing_lines);
dephasing_lines_temp = diff(dephasing_lines);
dephasing_diffs_xi = dephasing_lines_temp(:,:,[1,50,96]);

%% plotting
for d = 1:length(datadirs)
figure(1)
plot(dephasing_xi0_plot,squeeze(dephasing_diffs_all(:,d,:)) + dephasing_z,'linewidth',1.5);
xlabel('middle \xi_0 position for difference (cm)')
ylabel('difference (\xi_0(n+1) - \xi_0(n)) + z (\lambda/2)');
xlim([0 max(dephasing_xi0_plot)]);
title(['g = ',leg{d}]);
P.fig_handle = gcf;
P.plot_name = ['all_z',datadirs{d}];
% P.save_plot();

figure(2)
hold on
dephasing_diffs_xi_plot = squeeze(dephasing_diffs_xi(:,d,:));
for iz = 1:length([1,50,96])
    plot(dephasing_xi0_plot,dephasing_diffs_xi_plot(:,iz),'linewidth',2,...
    'color',colors{iz}); %,'linestyle',line_style{iz}
end
hold off

ylim([1.85,2.14]); %-----------------------------------------------

xlim([0 max(dephasing_xi0_plot)]);

xlabel('middle \xi_0 position for difference (cm)')
ylabel('difference (\xi_0(n+1) - \xi_0(n)) (\lambda/2)');
str_leg = 'z = ';
legend({[str_leg,num2str(dephasing_z(1),2),' m'],[str_leg,num2str(dephasing_z(50),2),' m'],...
    [str_leg,num2str(dephasing_z(96),2),' m']},'location','best')
title(['g = ',leg{d}]);
drawnow;
P.fig_handle = gcf;
P.plot_name = ['z15096',datadirs{d}];

stop = 1;
% P.save_plot();

figure(3)
plot(dephasing_z,squeeze(dephasing_lines(:,d,:)),'linewidth',1.5);
ylabel('zero - crossings (\lambda/2)')
xlabel('z (m)');
xlim([0 max(dephasing_z)]);
title(['g = ',leg{d}]);
P.fig_handle = gcf;
P.plot_name = ['all_z_deph',datadirs{d}];
P.save_plot();


close all
end % datadirs

%% plotting
% ind_xi = [53,25,1];
% for xi = 1:length(ind_xi)
%     
%     ax_x0(xi) = nexttile; 
%     ax_x0(xi).FontSize = fontsize_label;
%     
%     hold on
%     for d = 1:length(datadirs)
%         p1 = plot(dephasing_z,squeeze(dephasing_lines(ind_xi(xi),d,:)),...
%             line_style{d},'LineWidth',2);
%     end
% %     text(letterbox_x,letterbox_y,letters{i_letter},'FontSize',fontsize_label,'Units','normalized')
% %     i_letter = i_letter + 1;
%     
%     yline(-1,'--','LineWidth',1,'color',[0 0 0])
%     yline(-2,'--','LineWidth',1,'color',[0 0 0])
%     xline(4,'--','LineWidth',1,'color',[0 0 0]);
%     hold off
%     
%     switch dephasing_xi(ind_xi(xi))
%         case dephasing_xi(53)
%             position_word = '(back of the bunch)';
%         case dephasing_xi(25)
%             position_word = '(middle of the bunch)';
%         case dephasing_xi(1)
%             position_word = '(front of the bunch)';
%     end % dephasing xi (title)
%     title(['\xi_0 = ',num2str(dephasing_xi(ind_xi(xi))),' cm ',position_word]);
%     
%     xlim([0 10])
% %     ylim([-4 0.5]*2)
%     
% end % xi
% 
% legend(ax_x0(3),leg,'location','southeast','FontSize',fontsize_annotation,'NumColumns',2)
% ylabel(tt,'zero-crossing position shift (\lambda_p/2)','FontSize',fontsize_label)
% xlabel(tt,'z (m)','FontSize',fontsize_label)
% 
% P.fig_handle = fig_dephase;
% P.save_plot();
% 
% save('loading_files/gradsim_dephasing.mat','dephasing_z','dephasing_lines');



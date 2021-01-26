%________________________________________________________________________
% Script to produce calcualte the charge at some xi along z.
%
% Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 22/01/2020
%________________________________________________________________________

clear;
close all;

% load files
load('color_purple_to_green.mat');
% temp = load('sim_charge10.mat');
% sim_charge = temp.sim_charge;

% data directory
datadirs = {'gm10','g0','gp10'};
datadirs = {'gm20','gm15','gm10','gm5','g0','gp5','gp10','gp15','gp20'};
leg      = {'-1 %/m','0 %/m','1 %/m'};
leg = {'-2 %/m','-1.5 %/m','-1 %/m','-0.5 %/m','0 %/m','0.5 %/m','1 %/m','1.5 %/m','2 %/m'};

dataformat  = 'mat';
useAvg      = false;
initialdump = 0;
dump_list   = 0:5:100;

% save directory
plots_dir           = ['gradsim_test/fieldinxi/',''];
plot_name_suffix    = [''];
save_format         = {'png','eps','fig'};
plot_name           = 'fieldinxi';

% properties
plasma_density  = 1.81e14;
property        = 'fields';

trans_limit = 0.0536;
dephasing_xi = [14;7;1];


% Load the analysis class and initial charge
O = OsirisDenormalizer(...
    'datadir','g0','dataformat',dataformat,'useAvg',useAvg',...
    'dump',initialdump,'plasmaden',plasma_density,...
    'property',property,'wakefields_direction','trans',...
    'trans_range',[0,trans_limit]);
P = Plotty('plots_dir',plots_dir,'plasmaden',plasma_density,...
    'plot_name',plot_name,'save_flag',true);

xi_ranges = dephasing_xi + [O.plasma_wavelength/2,-O.plasma_wavelength/2];
fieldvsz = zeros(length(dephasing_xi),length(datadirs),length(dump_list));
plot_z = zeros(length(dephasing_xi),length(datadirs),length(dump_list));

% begin loop
for d = 1:length(datadirs)
    O.datadir = datadirs{d};
    
    
    for n = 1:length(dump_list)
        O.dump = dump_list(n); O.getdata(); O.assign_fields();
        
        for x = 1:length(dephasing_xi)
            
            O.xi_range = xi_ranges(x,:);
            
            %         if any(O.xi_range < 2) && (ismember(O.datadir,{'gm10','gp10'}))
            %             continue
            %         end
            
            O.denorm_Efield(); O.denorm_distance();
            
            z_ind = O.z > O.dtime+O.simulation_window - O.xi_range(1) & ... %large
                O.z <= O.dtime+O.simulation_window - O.xi_range(2); % small
            
            r_ind = O.r >= O.trans_range(1) & ...
                O.r < O.trans_range(2);
            trans_fieldinxi = O.transfield(r_ind,z_ind);
            i_field = trans_fieldinxi > 0;
            fieldvsz(x,d,n) = mean(trans_fieldinxi(i_field));
            plot_z(x,d,n) = O.propagation_distance/100;
            
        end % for xi
        O.progress_dump('dump',n,length(dump_list))
        
    end % for dump
    
    O.progress_dump('directory',d,length(datadirs))
end % for datadirs

i_color = [3,5,7];
i_color = [1:9];
%% plotting

fig_cvsz = figure(1);
fig_cvsz.OuterPosition = [100 100 1200 400];
line_style = {'--','-','-.'};
line_style = {'-','-','--','--','-','-.','-.',':',':'};

tt = tiledlayout(1,3);
tt.TileSpacing = 'compact';
tt.Padding = 'compact';
for x = 1:length(dephasing_xi)
    ax(x) = nexttile;
    for d = 1:length(datadirs)
        %         if (x == 3) && (ismember(datadirs{d},{'gm10','gp10'}))
        %             continue
        %         end
        hold on
        plot(squeeze(plot_z(x,d,:)),squeeze(fieldvsz(x,d,:)),...
            line_style{d},'LineWidth',2,'color',cc(i_color(d),:))
        hold off
    end % datadir
    hold on
    xline(4,'--','LineWidth',1,'color',[0 0.4470 0.7410]);
    hold off
    switch dephasing_xi(x)
        case 14
            position_word = '(back of the bunch)';
        case 7
            position_word = '(middle of the bunch)';
        case 1
            position_word = '(front of the bunch)';
    end
    title(['\xi_0 = ',num2str(dephasing_xi(x)),' cm ',position_word]);
    ylim([0,max(fieldvsz,[],'all')])
    xlim([0,10])
    
end % xi


legend(ax(3),leg,'location','northeast','FontSize',10)
xlabel(tt,'z (m)')
ylabel(tt,'mean defocusing field (MV/m)');


P.fig_handle = fig_cvsz;
P.save_plot();











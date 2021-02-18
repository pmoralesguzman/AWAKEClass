%________________________________________________________________________
% Plot the 2D wakefields together with the proton bunch, or each
% individually.
%
% Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 11/02/2020
%________________________________________________________________________

close all;
co =[    0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840
    0.9290    0.6940    0.1250];

% file location variables
datadirs = {'e+e-sqrt2','e+e-12','e+e-1sqrt2','e+e-1'};
legs = {'$\sigma_z = \sqrt{2}\,(1/k_{pe})$','$\sigma_z = 1/2\,(1/k_{pe})$',...
    '$\sigma_z = 1/\sqrt{2}\, (1/k_{pe})$','$\sigma_z = 1\,(1/k_{pe})$'};
linestyles = {'-',':','-.','--'};
dataformat = 'h5';
useAvg = 1;
dump_list = 2:1:2;

% saving data
save_flag = 0;
save_format = {'png','eps','fig'};

% plasma properties
plasmaden = 2e14;

% choose fields to plot
wakefields_direction = 'long'; % trans, long

% choose species density to plot
species = 'proton_beam';

% choose limits (in cm, must denormalize)
trans_range = [0 0.0536];
xi_range = [21 0];
lineout_point = 2;

% choose property to plot
property_plot = 'wakefields'; % density, wakefields, both

% create movie or not
create_movie = false;

% choose between normalized and denormalized units
denormalize_flag = false; % true, false

% choose if make pause or not
make_pause = false;

% directory to save the plots
plots_dir = ['positronplasma/','lineoutcompare/',...
    property_plot,'/',wakefields_direction,'/'];

P = Plotty('datadir',datadirs{1},'dataformat','h5',...
    'useAvg',useAvg,'dump_list',dump_list(1),...
    'save_flag',save_flag,'save_format',save_format,'plots_dir',plots_dir,...
    'create_movie',create_movie,...
    'plasmaden',plasmaden,'trans_range',trans_range,'xi_range',xi_range,...
    'wakefields_direction',wakefields_direction,'species',species,...
    'property_plot',property_plot,'denormalize_flag',denormalize_flag,...
    'make_pause',make_pause,'lineout_point',lineout_point);

for n = 1:length(dump_list)
    
    P.dump_list = dump_list(n);
    
    for d = 1:length(datadirs)
        P.datadir = datadirs{d};
        hold on
        P.lineout_plot();
        ax = gca;
        ax.ColorOrder = co;
        P.plot_handle.LineStyle = linestyles{d};
        hold off
    end
    
    legend(legs{:},'location','best','interpreter','latex')
%      legend(legs{:},'location','best')

    P.plot_name = ['lineout_comp','n',num2str(n)];
    P.fig_handle = gcf;
    ylabel('E_z (E_{WB})')
    title(['propagation dist. = ' num2str(P.propagation_distance,3),' cm']);
    P.save_flag = 1;
    P.save_plot();
end


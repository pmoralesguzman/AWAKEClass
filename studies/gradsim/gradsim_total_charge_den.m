%________________________________________________________________________
% Comparison between total charge measured experimentally and in
% simulations
%
% Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 06/10/2020
%________________________________________________________________________

clear;
% close all;
% load experimental data from Tatiana

blue = [47, 102, 169]/256;
red = [255, 64, 49]/256;

% charge_exp    = load('totalcharge_exp.txt');
% errorbars_exp = load('errorbars_totalcharge_exp.txt');

charge_exp    = load('totalcharge_1sigma.txt');
errorbars_exp = load('errorbars_totalcharge_1sigma.txt');

% -----------------------------------------------------------------
% data directory
datadirs    = {'gm20','gm15','gm10','gm5','g0','gp5','gp10','gp15','gp20'};
% datadirs    = {'gm20d2','gm15d2','gm10d2','gm5d2','g0d2','gp5d2','gp10d2','gp15d2','gp20d2'};

grads_exp   = [-19.4,-9.3,-5.16,0.3,4.3,8.7,13,20]/1;
grads_sim   = [-20,-15,-10,-5,0,5,10,15,20];
dataformat  = 'mat';
useAvg      = false;
initialdump = 0;
dump        = 119;

% save directory
plots_dir           = ['gradsim_paper/total_charge_compare/',''];
plot_name_suffix    = [''];
save_format         = {'png','eps','fig'};


% properties
plasma_density  = 1.81e14;

% plasma parameters
seeding_position    = 3.8074; % (127) ps ps*1e-12*O.c_cm
sigma_z             = 6.98516; % cm
sigma_exp           = 0.0536; % cm 0.0536, beta = 0.02*srqt(1+12^2/4.9^2) = 0.0529
exp_upperlimit      = sigma_exp;

% analysis
distance    = 200; % cm
limitr      = 1; %linspace(1,100,99)/100*2.9;

% switches
do_cvsr             = false; % flag to do individual c vs r plots
save_all_plots      = false;
save_plot_flag      = true;

% calculated variables
trans_limit = sigma_exp*limitr; % trans. limit in cm
exp_steps   = size(charge_exp,1);
exp_r       = linspace(sigma_exp/(exp_steps+1),3*exp_upperlimit,exp_steps);
% analytical charge fraction
bunch_trans_fraction= 1-exp(-(trans_limit/sigma_exp).^2/2);
bunch_long_fraction = normcdf(seeding_position,0,sigma_z);
bunch_back_fraction = normcdf(seeding_position - 3*sigma_z,0,sigma_z);
bunch_fraction_outside_simulation_window ...
    = bunch_trans_fraction*(1 - bunch_long_fraction + bunch_back_fraction);

% Load the analysis class and initial charge
O = OsirisDenormalizer('datadir','g0','dataformat',dataformat,'useAvg',useAvg',...
    'dump',initialdump,'plasmaden',plasma_density,...
    'property','density','raw_dataset','q',...
    'species','proton_beam','direction','r',...
    'trans_range',[0,trans_limit]);

plot_name = ['bunchfrac',num2str(dump),'s',num2str(limitr),'den'];
P = Plotty('plots_dir',plots_dir,'plasmaden',plasma_density,...
    'plot_name',plot_name,'save_flag',save_plot_flag);

% now that we have the analytical bunch fractions, calculate the initial charge
initial_charge      = 3e11;

% the charge within 1 sigma of the unmodulated proton bunch is taken as 1.0
% in the experiment
sigma_normalization = 1.0;
exp_normalization = 1 - exp(-(sigma_normalization)^2/2);

% initialize charges
charge_translim = zeros(length(trans_limit),length(grads_sim));
charge_fraction = zeros(length(trans_limit),length(grads_sim));


%% begin loop

O.dump = dump;

for d = 1:length(datadirs)
    
    % select study directory
    O.datadir = datadirs{d};
    O.getdata(); O.assign_density(); O.denorm_density(); O.denorm_distance();
    for r = 1:length(trans_limit)
        ind_r = O.r < trans_limit(r);        
        charge_translim(r,d) = O.cylindrical_integration(O.r(ind_r),O.z,O.proton_beam(ind_r,:),'sum');
        charge_fraction(r,d) = charge_translim(r,d)/initial_charge ...
            + bunch_fraction_outside_simulation_window(r);
        O.progress_dump('charge fraction r',r + (d-1)*length(trans_limit),length(trans_limit)*length(datadirs));
    end
    
end
sim_charge = charge_fraction/exp_normalization;

%% plot bunch population of selected bunch for the gradients
colors = {'r','k',[0 0.5 0],'b',[0.5,0,0.5]};
sz = 50;
marker_size = 10;

if length(trans_limit) == 1 % if only one radius is chosen, plot both on the same graph
    fig_cvsg = figure;
    [~,indr_exp] = min(abs(trans_limit-exp_r));
    plot_exp = errorbar(grads_exp/10,charge_exp(indr_exp,:),errorbars_exp(indr_exp,:),...
        's','MarkerSize',marker_size,'Color',blue); %    DATA
    set(plot_exp, 'MarkerFaceColor', get(plot_exp,'Color'))
    hold on
    plot_sim = plot(grads_sim/10,sim_charge,'o','MarkerSize',marker_size,'Color',red); % SIM
    set(plot_sim, 'MarkerFaceColor', get(plot_sim,'Color'))
%     p1 = plot(grads_exp/10,0.5595/93394892.27646959*[80703992.52017055 , 68925947.82344057 , 75240015.91677196 , 93394892.27646959 , 138102471.82438824 , 163936581.1924591 , 182302159.40055966 , 221130760.18843037],...
%         'hm');
%     p2 = plot(grads_exp/10,0.5595/10268414.871418467*[9179849.047179596, 8423956.737214955, 8853881.912626276, 10268414.871418467, 13686601.622400615, 15783961.96087067, 17334177.540318232, 20445446.223657124],...
%         'pk');
    
    hold off
%     legend([plot_sim(1) plot_exp(1) p1 p2],'Simulation','Experiment','SC charge','SC counts','location','northwest');
    legend([plot_sim(1) plot_exp(1)],'Simulation','Experiment','location','northwest');

    % 'Experiment \pm std. dev.'
    ax = gca;
    ax.FontSize = 14;
    grid on
    xlabel('density gradient (%/m)')
    ylabel('charge fraction')
    xlim([-2.2,2.2])
    %     ylim([0 0.85])
    
    
%     % experiment
%     x2 = [0.35,0.5];
%     y2 = [0.7,0.5];
%     a2 = annotation('textarrow',x2,y2,'LineWidth',1.5,'HeadStyle','none',...
%         'String','Experiment','FontSize',14);
%     a2.Color = blue;
%     
%     % simulation
%     x1 = [0.7,0.55];
%     y1 = [0.45,0.25];
%     a1 = annotation('textarrow',x1,y1,'LineWidth',1.5,'HeadStyle','none',...
%         'String','Simulation','FontSize',14);
%     a1.Color = red;
    
    drawnow;
    P.fig_handle = fig_cvsg;
    P.save_plot();
    
    do_cvsr = false;
    
else % otherwise, each set of curves in its own graph
    % plot of experimental data
    fig_fracallexp = figure(1);
    grads_exp_matrix = repmat(grads_exp/10,exp_steps,1);
    plot_exp = plot(grads_exp_matrix',charge_exp(1:exp_steps,:)','o-');
    xlabel('density gradient %/m')
    ylabel('charge fraction')
    title('charge fraction (experiment)')
    xlim([-2.2,2.2])
    ylim([0 1.2])
    drawnow;
    P.plot_name = ['allfracsexp',num2str(dump),''];
    P.fig_handle = fig_fracallexp;
    P.save_plot();
    
    % plot of simulation
    fig_fracallsim = figure(2);
    grads_mat = repmat(grads_sim/10,99,1);
    plot_sim = plot(grads_mat',sim_charge','-o');
    xlabel('density gradient %/m')
    ylabel('charge fraction')
    title('charge fraction (simulation)')
    xlim([-2.2,2.2])
    ylim([0 1.2])
    drawnow;
    P.plot_name = ['allfracsim',num2str(dump),''];
    P.fig_handle = fig_fracallsim;
    P.save_plot();
    
end % if length trans limit

if do_cvsr
    sim_charge(:,2) = [];
    grads_sim(2) = [];
    for g = 1:length(grads_exp)
        fig_cvsr = figure(g+100);
        hold on
        plot(trans_limit,charge_exp(:,g),'Linewidth',2)
        plot(trans_limit,sim_charge(:,g),'Linewidth',2)
        hold off
        xlabel('trans. limit')
        ylabel('charge fraction')
        title(['g = ',num2str(grads_sim(g)),' %'])
        legend('exp','sim','location','northwest')
        
        P.save_format = 'png';
        P.plot_name = ['cvsr',num2str(g)];
        P.fig_handle = fig_cvsr;
        P.save_plot();
    end
end



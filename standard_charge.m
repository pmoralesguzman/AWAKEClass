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

% clear;
close all;
% load experimental data from Tatiana

% charge_exp    = load('totalcharge_exp.txt');
% errorbars_exp = load('errorbars_totalcharge_exp.txt');

charge_exp    = load('totalcharge_1sigma.txt');
errorbars_exp = load('errorbars_totalcharge_1sigma.txt');

% data directory
datadirs    = {'gm20','gm15','gm10','gm5','g0','gp5','gp10','gp15','gp20'};
grads_exp   = [-19.4,-9.3,-5.16,0.3,4.3,8.7,13,20]/1;
grads_sim   = [-20,-15,-10,-5,0,5,10,15,20];
dataformat  = 'h5';
useAvg      = false;
initialdump = 0;
dump        = 100;

% save directory
plots_dir           = ['gradsim_paper/total_charge_compare/',''];
plot_name_suffix    = [''];
save_format         = {'png','eps','fig'};

% properties
plasma_density  = 1.81e14;

% bunch parameters
seeding_position    = 127; % ps
sigma_z             = 6.98516; % cm
sigma_exp           = 0.066; % cm 0.0536
exp_upperlimit      = sigma_exp;

% analysis
distance    = 350; % cm
limitr      = 1; %linspace(1,100,99)/100*3;

% switches
do_cvsr             = true; % flag to do individual c vs r plots
save_all_plots      = false;
save_plot_flag      = false;

% calculated variables
trans_limit = sigma_exp*limitr; % trans. limit in cm
exp_steps   = size(charge_exp,1);
exp_r       = linspace(sigma_exp/(exp_steps+1),3*exp_upperlimit,exp_steps);
% analytical charge fraction
seeding_position    = seeding_position*1e-12*O.c_cm; % cm
bunch_trans_fraction= 1-exp(-(trans_limit/sigma_exp).^2/2);
bunch_long_fraction = normcdf(seeding_position,0,sigma_z);
bunch_back_fraction = normcdf(seeding_position - 3*sigma_z,0,sigma_z);
bunch_fraction_outside_simulation_window ...
    = bunch_trans_fraction*(1 - bunch_long_fraction + bunch_back_fraction);

% Load the analysis class and initial charge

O = OsirisDenormalizer(...
    'datadir','g0','dataformat',dataformat,'useAvg',useAvg',...
    'dump',initialdump,'plasmaden',plasma_density,...
    'property','raw','raw_dataset','q',...
    'species','proton_beam','direction','r');
O.getdata(); O.assign_raw();

plot_name = ['bunchfrac',num2str(dump),'s',num2str(limitr),'raw'];
P = Plotty('plots_dir',plots_dir,'plasmaden',plasma_density,...
    'plot_name',plot_name,'save_flag',save_plot_flag,'save_format',save_format);

% now that we have the analytical bunch fractions, calculate the initial charge
initial_charge      = sum(O.q_raw)/(bunch_long_fraction - bunch_back_fraction);

% the charge within 1 sigma of the unmodulated proton bunch is taken as 1.0
% in the experiment
sigma_normalization = 1;
exp_normalization = 1 - exp(-(sigma_normalization)^2/2);
exp_normalization = 1*exp_normalization;

% initialize charges
charge_translim = zeros(length(trans_limit),length(grads_sim));
charge_fraction = zeros(length(trans_limit),length(grads_sim));


%% begin loop

O.dump = dump;

for d = 1:length(datadirs)
    
    % select study directory
    O.datadir = datadirs{d};
    
    O.raw_dataset = 'q'; O.getdata(); O.assign_raw();
    O.raw_dataset = 'x'; O.getdata(); O.assign_raw();
    O.raw_dataset = 'p'; O.getdata(); O.assign_raw();
    O.raw_dataset = 'ene'; O.getdata(); O.assign_raw();
    
    % push charges
    new_r = O.charge_pusher(O,distance);
    new_r = O.denorm_distance(new_r);
    
    for r = 1:length(trans_limit)
        ind_translim = new_r < trans_limit(r);
        charge_translim(r,d) = sum(O.q_raw(ind_translim))/initial_charge;
        charge_fraction(r,d) = charge_translim(r,d) ...
            + bunch_fraction_outside_simulation_window(r);
    end
    
    sim_charge = charge_fraction/exp_normalization;
    
end

%% plot bunch population of selected bunch for the gradients
colors = {'r','k',[0 0.5 0],'b',[0.5,0,0.5]};
sz = 50;
marker_size = 10;

if length(trans_limit) == 1 % if only one radius is chosen, plot both on the same graph
    fig_cvsg = figure(1);
    [~,indr_exp] = min(abs(trans_limit-exp_r));
    plot_exp = errorbar(grads_exp/10,charge_exp(indr_exp,:),errorbars_exp(indr_exp,:),...
        's','MarkerFaceColor','auto','MarkerSize',marker_size); %    DATA
    hold on
    plot_sim = plot(grads_sim/10,sim_charge,'o','MarkerSize',marker_size); % SIM
    set(plot_sim, 'MarkerFaceColor', get(plot_sim,'Color'))
    hold off
    legend([plot_sim(1) plot_exp(1)],'Simulation','Exp. data \pm std. dev.','location','northwest');
    grid on
    xlabel('density gradient (%/m)')
    ylabel('total charge fraction (a.u.)')
    xlim([-2.2,2.2])
    %     ylim([0 0.85])
    drawnow;
    P.fig_handle = fig_cvsg;
    P.save_plot();
    
    do_cvsr = false;
    
else % otherwise, each set of curves in its own graph
    % plot of experimental data
    fig_fracallexp = figure(1);
    grads_exp_matrix = repmat(grads_exp/10,exp_steps,1);
    plot_exp = plot(grads_exp_matrix',1/1*charge_exp(1:exp_steps,:)','o-');
    xlabel('density gradient %/m')
    ylabel('bunch fraction')
    title('Total charge (experiment)')
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
    ylabel('bunch fraction')
    title('Total charge (simulation)')
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



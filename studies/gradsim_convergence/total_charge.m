%________________________________________________________________________
% Comparison between total for simulations with different input parameters
% to test convergence
%
% Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 06/09/2030
%________________________________________________________________________

% clear;
close all;

%% input info

% directories where the files are
datadirs = {'g0','g0zh','g0rh','g0dt92','g0dt95'};

% plasma parameters
plasma_density = 1.81e14;
x_ticks = [1,2,3];
seeding_position = 127; % ps
sigma_z = 6.98516; % cm
sigma_exp = 0.066; % cm 0.0536
exp_upperlimit = sigma_exp;

% directory to save plots
plots_dir   = 'total_charge_convergence';

% property to be analized
property    = 'density';
species     = 'proton_beam';
do_cvsr     = true; % flag to do individual c vs r plots

% simulation parameters
initialdump = 0;
studydump = 100;
limitr = 1;
dataformat = 'h5';

% calculated variables
trans_limit = sigma_exp*limitr; % trans. limit in cm
exp_steps   = size(charge_exp,1);
exp_r       = linspace(sigma_exp/(exp_steps+1),3*exp_upperlimit,exp_steps);

% Load the analysis class and initial charge
O = OsirisDenormalizer('plasmaden',1.81e14,'property','density',...
    'dump',initialdump,'dataformat',dataformat,'species','proton_beam',...
    'trans_range',[0,trans_limit]);

% analytical charge fraction
seeding_position    = seeding_position*1e-12*O.c_cm; % cm
bunch_trans_fraction= 1-exp(-(trans_limit/sigma_exp).^2/2);
bunch_long_fraction = normcdf(seeding_position,0,sigma_z);
bunch_back_fraction = 0*normcdf(seeding_position - 3*sigma_z,0,sigma_z);
bunch_fraction_outside_simulation_window ...
                    = bunch_trans_fraction*(1 - bunch_long_fraction + bunch_back_fraction);

% now that we have the analytical bunch fractions, calculate the initial charge
initial_charge      = 3e11;

% the charge within 1 sigma of the unmodulated proton bunch is taken as 1.0
% in the experiment
sigma_normalization = 0.94;  
exp_normalization = 1 - exp(-(sigma_normalization)^2/2);
% exp_normalization = 1;

% initialize charges
charge_translim = zeros(length(trans_limit),length(x_ticks));
charge_fraction = zeros(length(trans_limit),length(x_ticks));


%% begin loop

O.dump = studydump;

for d = 1:length(datadirs)
    
    % select study directory
    O.datadir = datadirs{d};
    O.getdata(); O.trim_data(); O.denorm_density();
    
    for r = 1:length(trans_limit)
        charge_translim(r,d) = O.cylindrical_integration(O.r,O.z,O.proton_beam,'simpsons')/initial_charge;
        charge_fraction(r,d) = charge_translim(r,d) ...
            + bunch_fraction_outside_simulation_window(r);
    end
    
    sim_charge = charge_fraction/exp_normalization;
    
end

%% plot bunch population of selected bunch for the gradients
% figure(1)
colors = {'r','k',[0 0.5 0],'b',[0.5,0,0.5]};
sz = 50;

if length(trans_limit) == 1 % if only one radius is chosen, plot both on the same graph
    
    plot_sim = scatter(x_ticks,sim_charge,sz,'o','filled');
    grid on
    xlabel('simulation')
    ylabel('bunch fraction')
    %     ylim([0 0.85])
    set(gca,'xtick',x_ticks,'xticklabel',names)
    drawnow;
    P = Plotty('fig_handle',gcf,'plots_dir',plots_dir,'plot_name',['bunchfrac',num2str(studydump),'']);
    P.save_plot();
    
    do_cvsr     = false;
    
else % otherwise, each set of curves in its own graph
    
    
    % plot of simulation
    figure(2)
    grads_mat = repmat(x_ticks,99,1);
    plot_sim = plot(grads_mat',sim_charge','-o');
    xlabel('simulation')
    ylabel('bunch fraction')
    ylim([0 1.2])
    
    drawnow;
    P = Plotty('fig_handle',gcf,'plots_dir',plots_dir,'plot_name',['allfracsim',num2str(studydump),'']);
    P.save_plot();
    
end % if length trans limit

if do_cvsr
    for d = 2:length(datadirs)
        figure(g+100)
        hold on
        plot(trans_limit,sim_charge(:,g),'Linewidth',2)
        plot(trans_limit,sim_charge(:,1),'Linewidth',2)
        hold off
        
        xlabel('trans. limit')
        ylabel('charge fraction')
        title(['','Charge transverse cumulative distribution'])
        legend('g0',datadirs{d},'location','northwest')
        P = Plotty('plots_dir',plots_dir,'plot_name',['cvsr',datadirs{d}],...
            'save_format','png',fig_handle,gcf);
        P.save_plot();
    end
end



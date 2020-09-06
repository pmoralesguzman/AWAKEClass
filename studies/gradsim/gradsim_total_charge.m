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
% Last update: 04/06/2030
%________________________________________________________________________

% clear;
close all;

%% input info

% load experimental data from Tatiana
charge_exp    = load('totalcharge_1sigma.txt');
errorbars_exp = load('errorbars_totalcharge_1sigma.txt');

% charge_exp    = load('totalcharge_exp.txt');
% errorbars_exp = load('errorbars_totalcharge_exp.txt');


% directories where the files are
datadirs = {'gm20','gm15','gm10','gm5','g0','gp2','gp5','gp10','gp15','gp20'};
% datadirs = {'real_studies/gm20'};

% plasma parameters
plasma_density = 1.81e14;
grads_exp = [-19.4,-9.3,-5.16,0.3,4.3,8.7,13,20]/1;
grads_sim = [-20,-15,-10,-5,0,2,5,10,15,20];
seeding_position = 127; % ps
sigma_z = 6.98516; % cm
sigma_exp = 0.0536; % cm 0.0536
exp_upperlimit = sigma_exp;

% directory to save plots
plots_dir   = 'total_charge_compare';

% property to be analized
property    = 'density';
species     = 'proton_beam';
do_cvsr     = true; % flag to do individual c vs r plots

% simulation parameters
initialdump = 0;
studydump = 100;
limitr = 1; %linspace(1,100,99)/100*3;
dataformat = 'h5';
distance = 350; % cm

% calculated variables
trans_limit = sigma_exp*limitr; % trans. limit in cm
exp_steps   = size(charge_exp,1);
exp_r       = linspace(sigma_exp/(exp_steps+1),3*exp_upperlimit,exp_steps);

% Load the analysis class and initial charge
O = OsirisDenormalizer('plasmaden',1.81e14,'property','raw','direction','r',...
    'dump',initialdump,'dataformat',dataformat,'species','proton_beam');
O.datadir = 'g0'; O.raw_dataset = 'q'; O.getdata(); O.assign_raw();

% analytical charge fraction
seeding_position    = seeding_position*1e-12*O.c_cm; % cm
bunch_trans_fraction= 1-exp(-(trans_limit/sigma_exp).^2/2);
bunch_long_fraction = normcdf(seeding_position,0,sigma_z);
bunch_back_fraction = normcdf(seeding_position - 3*sigma_z,0,sigma_z);
bunch_fraction_outside_simulation_window ...
                    = bunch_trans_fraction*(1 - bunch_long_fraction + bunch_back_fraction);

% now that we have the analytical bunch fractions, calculate the initial charge
initial_charge      = sum(O.q_raw)/(bunch_long_fraction - bunch_back_fraction);

% the charge within 1 sigma of the unmodulated proton bunch is taken as 1.0
% in the experiment
sigma_normalization = 1;  
exp_normalization = 1 - exp(-(sigma_normalization)^2/2);
exp_normalization = 0.85*exp_normalization;

% initialize charges
charge_translim = zeros(length(trans_limit),length(grads_sim));
charge_fraction = zeros(length(trans_limit),length(grads_sim));


%% begin loop

O.dump = studydump;

for d = 1:length(datadirs)
    
    % select study directory
    datadir = datadirs{d};
    O.datadir = datadir;
    
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
% figure(1)
colors = {'r','k',[0 0.5 0],'b',[0.5,0,0.5]};
sz = 50;

if length(trans_limit) == 1 % if only one radius is chosen, plot both on the same graph
    [~,indr_exp] = min(abs(trans_limit-exp_r));
    plot_exp = errorbar(grads_exp/10,charge_exp(indr_exp,:),errorbars_exp(indr_exp,:),...
        'o','MarkerFaceColor','auto');
    hold on
    plot_sim = scatter(grads_sim/10,sim_charge,sz,'o','filled');
    hold off
    legend([plot_sim(1) plot_exp(1)],'Simulation','Exp. data','location','northwest');
    grid on
    xlabel('density gradient %/m')
    ylabel('bunch fraction')
    xlim([-2.2,2.2])
%     ylim([0 0.85])
    drawnow;
    P = Plotty('fig_handle',gcf,'plots_dir',plots_dir,'plot_name',['bunchfrac',num2str(studydump),'both']);
   P.save_plot();
    
    do_cvsr     = false; 
    
else % otherwise, each set of curves in its own graph
    % plot of experimental data
    grads_exp_matrix = repmat(grads_exp/10,exp_steps,1);
    plot_exp = plot(grads_exp_matrix',1/1*charge_exp(1:exp_steps,:)','o-');
    xlabel('density gradient %/m')
    ylabel('bunch fraction')
    title('Total charge (experiment)')
    xlim([-2.2,2.2])
    ylim([0 1.2])
    drawnow;
    P = Plotty('fig_handle',gcf,'plots_dir',plots_dir,'plot_name',['allfracsexp',num2str(studydump),'']);
    P.save_plot();
    
    % plot of simulation
    figure(2)
    grads_mat = repmat(grads_sim/10,99,1);
    plot_sim = plot(grads_mat',sim_charge','-o');
    xlabel('density gradient %/m')
    ylabel('bunch fraction')
    title('Total charge (simulation)')
    xlim([-2.2,2.2])
    ylim([0 1.2])
    drawnow;
    P = Plotty('fig_handle',gcf,'plots_dir',plots_dir,'plot_name',['allfracsim',num2str(studydump),'']);
    P.save_plot();
    
end % if length trans limit

if do_cvsr
    sim_charge(:,2) = [];
    grads_sim(2) = [];
    for g = 1:length(grads_exp)
        figure(g+100)
        hold on
        plot(trans_limit,charge_exp(:,g),'Linewidth',2)
        plot(trans_limit,sim_charge(:,g),'Linewidth',2)
        hold off
        xlabel('trans. limit')
        ylabel('charge fraction')
        title(['g = ',num2str(grads_sim(g)),' %'])
        legend('exp','sim','location','northwest')
        P = Plotty('plots_dir',plots_dir,'plot_name',['cvsr',num2str(g)],'save_format','png');
        P.fig_handle = gcf;
        P.save_plot();
    end
end


%% functions in script
% function new_r = charge_pusher(OD,distance)
% 
% norm_time = OD.norm_distance(distance);
% 
% vel = OD.npr_raw./(OD.nE_raw+1);
% new_r = abs(OD.nr_raw + vel.*norm_time);
% 
% end


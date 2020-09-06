%________________________________________________________________________
% Evolution of the total charge along propagation for each gradient
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

clear;
close all;

%% input info

% directories where the files are
% datadirs = {'gm20','gm15','gm10','gm5','g0','gp2','gp5','gp10','gp15','gp20'};
datadirs = {'gm20','gm15','gm10','gm5','g0','gp5','gp10','gp15','gp20'};

% datadirs = {'gm20'};
load('color_purple_to_green.mat');
temp = load('sim_charge10.mat');
sim_charge = temp.sim_charge;

% plasma parameters
plasma_density = 1.81e14;
grads_sim = [-20,-15,-10,-5,0,5,10,15,20];
seeding_position = 127; % ps
sigma_z = 6.98516; % cm
% sigma_exp = 0.0536; % cm
sigma_exp = 0.066;
exp_upperlimit = sigma_exp;

% directory to save plots
plots_dir   = 'charge_evolution';
dataformat = 'mat';
useAvg = false;

% property to be analized
property    = 'density';
species     = 'proton_beam';

% simulation parameters
dump_initial_charge = 0;
dump_list = 0:1:134;
limitr = 1;%/0.0536;

distance = 350; % cm

% calculated variables
trans_limit = sigma_exp*limitr; % trans. limit in cm

% Load the analysis class and initial charge
O = OsirisDenormalizer('datadir','gm20','plasmaden',1.81e14,'property',property,'direction','r',...
    'dump',dump_initial_charge,'dataformat',dataformat,'species','proton_beam',...
    'raw_dataset','q','useAvg',useAvg,'trans_range',[0,trans_limit]);
O.getdata(); O.denorm_distance(); O.assign_density();  O.denorm_density();

% analytical charge fraction
seeding_position    = seeding_position*1e-12*O.c_cm; % cm
bunch_trans_fraction= (normcdf(trans_limit,0,sigma_exp)-normcdf(-trans_limit,0,sigma_exp))^2;
bunch_long_fraction = normcdf(seeding_position,0,sigma_z);
bunch_back_fraction = 0*normcdf(seeding_position - 3*sigma_z,0,sigma_z);
bunch_fraction_outside_simulation_window ...
    = (bunch_trans_fraction*(1 - bunch_long_fraction) + bunch_back_fraction);

% now that we have the analytical bunch fractions, calculate the initial charge
% initial_charge      = sum(O.q_raw)/(bunch_long_fraction - bunch_back_fraction);
initial_charge = O.cylindrical_integration(O.r,O.z,O.(O.species),'simpsons')/(bunch_long_fraction - bunch_back_fraction);

% the charge within 1 sigma of the unmodulated proton bunch is taken as 1.0
% in the experiment
sigma_normalization = 0.95;
exp_normalization = 1 - exp(-(sigma_normalization)^2/2);
% exp_normalization = sigma_normalization*exp_normalization;

% initialize charges
charge_translim = zeros(length(trans_limit),length(grads_sim));
charge_fraction = zeros(length(trans_limit),length(grads_sim));


%% begin loop
current_charge = zeros(length(datadirs),length(dump_list));
z = zeros(length(datadirs),length(dump_list));
for d = 1:length(datadirs)
    
    % select study directory
    datadir = datadirs{d};
    O.datadir = datadir;
    
    
    for n = 1:length(dump_list)
        O.dump = dump_list(n); O.getdata(); 
        O.trim_data(); O.denorm_density();
        switch O.property
            case 'raw'
                O.raw_dataset = 'q'; O.getdata(); O.assign_raw();
                O.raw_dataset = 'x'; O.getdata(); O.assign_raw();
                O.raw_dataset = 'p'; O.getdata(); O.assign_raw();
                O.raw_dataset = 'ene'; O.getdata(); O.assign_raw();
                
                % push charges
                old_r = O.nr_raw;
                new_r = charge_pusher(O,distance);
                new_r = O.denorm_distance(new_r);
                
                for r = 1:length(trans_limit)
                    ind_translim = new_r < trans_limit(r);
                    charge_translim(r,d) = sum(O.q_raw(ind_translim))/initial_charge;
                    charge_fraction(r,d) = charge_translim(r,d) ...
                        + bunch_fraction_outside_simulation_window(r);
                end
            case 'density'
                current_charge(d,n) = O.cylindrical_integration(O.r,O.z,O.(O.species));
                z(d,n) = O.propagation_distance();
                O.progress_dump(['charge evo. ',datadir],n,length(dump_list));
        end
    end
    
end

charge_fraction0 = (current_charge/initial_charge + ...
    bunch_fraction_outside_simulation_window)/exp_normalization;
% charge_fraction = [charge_fraction0,sim_charge'];
charge_fraction = charge_fraction0;
z = [z/100];

%% plotting 

figc = figure(1);
colororder(cc);
plot(z',charge_fraction','LineWidth',2)
xlim([0,13.5])

xlabel('propagation direction (m)')
ylabel('charge fraction (a.u.)');
legend('gm20','gm15','gm10','gm5','g0','gp5','gp10','gp15','gp20','Location','best')
% legend('gm20','gm15','gm10','gm5','g0','gp2','gp5','gp10','gp15','gp20','Location','best')
P = Plotty('plots_dir','charge_evolution','plot_name','gradsim_135','fig_handle',figc);
P.save_plot();


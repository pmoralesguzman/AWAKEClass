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

answer = questdlg('Are you sure?');

switch answer
    case 'Yes'
        % nothing
    otherwise
        return;
end


% clear;
close all;

run_it = 0;

% load files
load('color_purple_to_green.mat');
load('color_red_to_blue.mat');

% temp = load('sim_charge10.mat');
% sim_charge = temp.sim_charge;

% data directory
datadirs    = {'gm20','gm15','gm10','gm5','g0','gp5','gp10','gp15','gp20'};
% datadirs    = {'gm20'};

grads_sim   = [-20,-15,-10,-5,0,5,10,15,20];
leg         = {'-2 %/m','-1.5 %/m','-1 %/m','-0.5 %/m','0 %/m','0.5 %/m','1 %/m','1.5 %/m','2 %/m'};
dataformat  = 'mat';
useAvg      = false;
initialdump = 0;
dump_list   = 0:1:120;

% save directory
plots_dir           = ['gradsim_paper/charge_evolution/',''];
plot_name_suffix    = [''];
save_format         = {'png','eps','fig'};
plot_name           = 'chargeevolutionn120';

% properties
plasma_density  = 1.81e14;
property        = 'density';

% plasma parameters
seeding_position    = 3.8074; % (127) ps ps*1e-12*O.c_cm
sigma_z             = 6.98516; % cm
% sigma_exp = 0.0536; % cm
sigma_exp           = 0.0536; % cm 0.0536
exp_upperlimit      = sigma_exp;

% analysis
limitr      = 1; %linspace(1,100,99)/100*3;

% switches
save_plot_flag      = true;

% calculated variables
trans_limit = sigma_exp*limitr; % trans. limit in cm


% analytical charge fraction

bunch_long_fraction = normcdf(seeding_position,0,sigma_z);
bunch_back_fraction = 0*normcdf(seeding_position - 3*sigma_z,0,sigma_z);


% now that we have the analytical bunch fractions, calculate the initial charge
initial_charge = 3e11;



% initialize charges
chargevsz = zeros(length(datadirs),length(dump_list));
plot_z = zeros(length(datadirs),length(dump_list));
charge_translim = zeros(length(datadirs),length(dump_list));
charge_fraction = zeros(length(datadirs),length(dump_list));
exp_normalization = zeros(length(datadirs),length(dump_list));
sigma_beam = zeros(length(datadirs),length(dump_list));
bunch_trans_fraction = zeros(length(datadirs),length(dump_list));
bunch_fraction_outside_simulation_window = zeros(length(datadirs),length(dump_list));

% Load the analysis class and initial charge
O = OsirisDenormalizer(...
    'datadir','g0','dataformat',dataformat,'useAvg',useAvg',...
    'dump',initialdump,'plasmaden',plasma_density,...
    'property',property,'raw_dataset','q',...
    'species','proton_beam','direction','r',...
    'trans_range',[0,trans_limit]);

P = Plotty('plots_dir',plots_dir,'plasmaden',plasma_density,...
    'plot_name',plot_name,'save_flag',save_plot_flag);

if run_it ==1
    
    % begin loop
    for d = 1:length(datadirs)
        
        % select study directory
        
        O.datadir = datadirs{d};
        
        for n = 1:length(dump_list)
            O.dump = dump_list(n); O.getdata();
            
            switch O.property
                case 'raw'
                    
                    O.raw_dataset = 'x'; O.getdata(); O.assign_raw();
                    O.r_raw = O.denorm_distance(O.nr_raw);
                    O.raw_dataset = 'q'; O.getdata(); O.assign_raw();
                    O.denorm_distance();
                    
                    ind_translim = O.r_raw < trans_limit;
                    if n == 1
                        if dump_list(n) < 2
                            initial_charge = sum(O.q_raw)/(bunch_long_fraction - bunch_back_fraction);
                        else
                            error('for raw the first dump must be 0 or 1 to calculate initial charge')
                        end % if dump list 2
                    end % if n 1
                    
                    chargevsz(d,n) = sum(O.q_raw(ind_translim));
                    O.propagation_distance = O.dtime - 21/2;
                    
                case 'density'
                    O.trim_data(); O.denorm_density();
                    chargevsz(d,n) = O.cylindrical_integration(O.r,O.z,O.(O.species),'trapz');
                    
            end %switch property
            plot_z(d,n) = O.propagation_distance/100;
            
            % the charge within 1 sigma of the unmodulated proton bunch is taken as 1.0
            % in the experiment
            sigma_beam(d,n) = 0.02*sqrt(1+plot_z(d,n)^2/4.9^2);
            sigma_normalization = sigma_exp/sigma_beam(d,n);
            exp_normalization(d,n) = 1 - exp(-(sigma_normalization)^2/2);
            
            bunch_trans_fraction(d,n) = 1-exp(-(trans_limit./sigma_beam(d,n))^2/2);
            bunch_fraction_outside_simulation_window(d,n) = ...
                bunch_trans_fraction(d,n)*(1 - bunch_long_fraction + bunch_back_fraction);
            
            O.progress_dump(['charge evo. ',O.datadir],n,length(dump_list));
        end % for dump list
        
    end % for datadirs
    
    
    
    
    charge_fraction0 = (chargevsz/initial_charge + ...
        bunch_fraction_outside_simulation_window)./exp_normalization;
    % charge_fraction = [charge_fraction0,sim_charge'];
    charge_fraction = charge_fraction0;
    plot_z = [plot_z];
    
else
    
    load('loading_files/gradsim_charge_evolution.mat');
    
end %if run it

%% plotting

fig_cvsz = figure(1);
colororder(ccrb);
% line_style = {'-.','-.','-.','-.','-','--','--','--','--'};
line_style = {':','--','-.','-','-','-','-.','--',':'};
fontsize_annotation = 12;
fontsize_label = 14;

for d = 1:length(datadirs)
    hold on
    plot(plot_z(d,:),charge_fraction(d,:),line_style{d},'LineWidth',2)
    hold off
end
hold on
xline(4,'--','LineWidth',1,'color',[0 0.4470 0.7410]);
xline(10,'--','LineWidth',1,'color',[0 0.0 0.0]);
hold off
xlim([0,12])
ax = gca;
ax.FontSize = fontsize_label;
xlabel('z (m)','FontSize',fontsize_label)
ylabel('charge fraction','FontSize',fontsize_label);
legend(leg,'Location','southwest','Autoupdate','off','FontSize',fontsize_annotation)

% Create textarrow
annotation(fig_cvsz,'textarrow',[0.621230158730159 0.576785714285714],...
    [0.412169312169312 0.418518518518519],'String','g = 0 \%/m',...
    'Interpreter','latex',...
    'HeadStyle','none','FontSize',fontsize_annotation);

% Create textarrow
annotation(fig_cvsz,'textarrow',[0.574404761904762 0.606150793650794],...
    [0.794179894179895 0.526455026455027],'String','positive gradient values',...
    'Interpreter','latex',...
    'Color',ccrb(8,:),...
    'HeadStyle','none','FontSize',fontsize_annotation);

% Create textarrow
annotation(fig_cvsz,'textarrow',[0.537722222222222 0.481053968253968],...
    [0.17384126984127 0.354851851851852],'String','negative gradient values',...
    'Interpreter','latex',...
    'Color',ccrb(2,:),...
    'HeadStyle','none','FontSize',fontsize_annotation);

drawnow;

P.fig_handle = fig_cvsz;
P.save_plot();

save('loading_files/gradsim_charge_evolution.mat','plot_z','charge_fraction');

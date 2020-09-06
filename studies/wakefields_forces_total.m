%__________________________________________________________________________
% Script that calculates the amount of wakefields forces acting on
% particles (wakefields field * charge)
%
% For use with: Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 04/06/2020
%__________________________________________________________________________
clear;

% data directory
datadirs = {'gm20','gm15','gm10','gm5','g0','gp5','gp10','gp15','gp20'};
% datadirs = {'gm20'};

grads_sim = [-20,-15,-10,-5,0,5,10,15,20];

% parameters
plasma_density = 1.81e14;

% properties
property = 'fields';
species = 'proton_beam';
field = 'e';
direction = 'r';
wakefields_firection = 'trans';

% simulation parameters
dump_list = 1:1:100;
dataformat = 'h5';
useAvg = true;

OD = OsirisDenormalizer('datadir','gm20',...
    'plasmaden',plasma_density,'property',property,'species',species,...
    'field',field,'direction',direction,...
    'trans_range',[0.00 0.15],...
    'wakefields_direction',wakefields_firection,...
    'dump',0,'dataformat',dataformat,'useAvg',useAvg);

defocusing_force_total = zeros(1,length(datadirs));
focusing_force_total = zeros(1,length(datadirs));

for d = 1:length(datadirs)
    datadir = datadirs{d};
    
    for n = 1:length(dump_list)
        OD.datadir = datadirs{d};
        OD.dump = dump_list(n); OD.property = 'fields'; 
%         
        OD.getdata(); OD.denorm_distance(); OD.trim_data();
        OD.denorm_Efield(); 
        
        % OD.transfield, OD.r, OD.z
        
        OD.property = 'density'; OD.getdata(); OD.denorm_distance(); OD.trim_data();
        OD.denorm_density(); % OD.proton_beam
        proton_charge = (2*pi*OD.dz*diff([OD.r,OD.r(end)+OD.dr].^2)').*OD.proton_beam;
%         OD.property = 'raw';
%         OD.raw_dataset = 'q'; OD.getdata(); OD.assign_raw(); % OD.q_raw
%         
%         OD.raw_dataset = 'x';
%         OD.direction = 'r';
%         OD.getdata(); OD.assign_raw(); % OD.nr_raw
%         OD.r_raw = OD.denorm_distance(OD.nr_raw);
%         
%         OD.direction = 'z';
%         OD.getdata(); OD.assign_raw(); % OD.nr_raw
%         OD.z_raw = OD.denorm_distance(OD.nz_raw);
%         
%         % bin charge
%         
%         OD.dz;
%         
%         z_ind = ceil((OD.z_raw-OD.dtime)./OD.dz);
%         r_ind = ceil((OD.r_raw)./OD.dr);
%         
%         for ii = 1:length(z_ind)
%         Q(r_ind(ii),z_ind(ii)) = Q(r_ind(ii),z_ind(ii)) + OD.q_raw(ii);
%         end

        trans_force = OD.transfield.*proton_charge;
        defocusing_force = sum(trans_force(trans_force>0));
        trans_force = trans_force(OD.r>0.02,:);
        focusing_force = sum(trans_force(trans_force<0));
%         neutral_force =  sum(trans_force(trans_force>0));
        defocusing_force_total(d) = defocusing_force_total(d) + defocusing_force;
        focusing_force_total(d) = focusing_force_total(d) + focusing_force;
        
        
        
        fprintf('%.d / %.d ; %.d / %.d \n',n,length(dump_list),d,length(datadirs))
    end % end for dump list
    
    
end % for datadirs

fig1 = figure(1);
bar(grads_sim,defocusing_force_total);
xlabel('gradient');
ylabel('total defocusing force (a. u.)')

P = Plotty('fig_handle',fig1,'plot_name','bar_def','plots_dir','force');
P.save_plot();

fig2 = figure(2);
bar(grads_sim,focusing_force_total);
xlabel('gradient');
ylabel('total focusing force (a. u.)')

P = Plotty('fig_handle',fig2,'plot_name','bar_foc','plots_dir','force');
P.save_plot();

fig3 = figure(3);
diff_total = focusing_force_total + defocusing_force_total;
bar(grads_sim,diff_total);
xlabel('gradient');
ylabel('total focusing - defocusing force (a. u.)')

P = Plotty('fig_handle',fig3,'plot_name','bar_diff','plots_dir','force');
P.save_plot();

save('save_files/forces.mat','focusing_force_total','defocusing_force_total',...
    'diff_total');









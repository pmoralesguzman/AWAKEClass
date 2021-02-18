%________________________________________________________________________
% Create a time density profile of the proton bunch, similar to what is obtained in
% the experiment (density profile at one point in space, while time
% advances, instead of a profile frozen in time while z varies).
%
% Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 30/11/2020
%________________________________________________________________________

clear;
close all;

% maximum delay in the experiment: 618.2097 ps --> 18.5335 cm
% longitudinal pixels after the seeding: 1500
% to get 1500 pixel res., binsize = 0.41214 ps

% exp. total transverse distance = 672*0.00217 = 1.4582 cm
% exp. half distance ~= 0.7291 cm in ~336 pixels
% --> exp pixels in 0.16 cm -> 
% transverse pixels total = 672 --> 73.73 ~ 74

% parameters

datadir = 'gm10';
measurement_point = 1200; % in cm, from plasma beginning
dump = 0;
binsize = 0.41214;
trans_lims = linspace(0,0.6,2*138); %cm (74 for 0.16 cm (see above))
% trans_lims = linspace(0,0.01,40); %cm (74 for 0.16 cm (see above))

xi_range = [18.5335,0];
dataformat = 'h5';

% get macroparticle info

O = OsirisDenormalizer('plasmaden',1.81e14,'datadir',datadir,...
    'property','raw','xi_range',xi_range,...
    'dump',dump,'dataformat',dataformat,'species','proton_beam');
P = Plotty('plasmaden',1.81e14,'datadir',datadir,...
    'dump',dump,'dataformat',dataformat,'species','proton_beam',...
    'property_plot','density');
O.raw_dataset = 'x'; O.direction = 'z'; O.getdata(); O.assign_raw();
O.raw_dataset = 'x'; O.direction = 'r'; O.getdata(); O.assign_raw();


O.raw_dataset = 'p'; O.direction = 'z'; O.getdata(); O.assign_raw();
% pz = O.denorm_distance(O.npz_raw);
O.raw_dataset = 'p'; O.direction = 'r'; O.getdata(); O.assign_raw();
% pr = O.denorm_distance(O.npr_raw);

O.raw_dataset = 'ene'; O.getdata(); O.assign_raw();
O.raw_dataset = 'q'; O.getdata(); O.assign_raw();


% to get simulation window size, needed for the distance
O.dataformat = 'mat'; % done for gm20
O.property = 'density'; O.getdata(); O.trim_data(); O.denorm_distance(); 
O.dataformat = 'h5';

% manual trim of raw data

ixi = (O.dtime + O.simulation_window - O.denorm_distance(O.nz_raw)) < xi_range(1);

z = O.denorm_distance(O.nz_raw(ixi));
r = O.denorm_distance(O.nr_raw(ixi));
O.npr_raw = O.npr_raw(ixi); 
O.nr_raw = O.nr_raw(ixi);
O.q_raw = O.q_raw(ixi);
q = O.q_raw;
O.nE_raw = O.nE_raw(ixi);

% front and back propagate

prop_distances = measurement_point + O.simulation_window - z;

new_r = O.denorm_distance(O.charge_pusher(O,prop_distances));

delay_time = (O.dtime + O.simulation_window - z)/O.c_cm*1e12; % ps
t_simulation_window = (O.simulation_window)/O.c_cm*1e12; % ps

trim_window_size = (max(O.z)-min(O.z))/O.c_cm*1e12;

chargematrix = zeros(length(trans_lims)-1,ceil(trim_window_size/binsize));

for rr = 1:(length(trans_lims)-1)
    ir = (new_r >= trans_lims(rr)) & (new_r < trans_lims(rr+1));
    if sum(ir) == 0; continue; end
        
    
    
    % select only those particles inside the transverse limits
    q_r = O.q_raw(ir);
    delay_time_r = delay_time(ir);
    [N,edges,ind_bin] = histcounts(delay_time_r,0:binsize:trim_window_size);
    A = accumarray(ind_bin+1,delay_time_r);
    % bin
    ind_bin = ceil(delay_time_r/binsize);
    
    % sort the indeces of the binning and then sort the
    % elements in q
    [ind_sort,ind_order] = sort(ind_bin);
    q_sort = q_r(ind_order);
    
    % say how many there is in each bin in a cumulative way,
    % to use that as indices
    % for q, to quickly build up the indexed density along z
    ind_sum = [0,0,histcounts(ind_sort,1:1:max(ind_bin),'Normalization','cumcount')];
    
    % go through each element in bz and arrange that charge
    % according to the bins, using the indeces from histconts
    for bz = 1:max(ind_bin)
        chargematrix(rr,bz) = sum(q_sort(ind_sum(bz)+1:ind_sum(bz+1)));
    end
    O.progress_dump('building density',rr,length(trans_lims)-1)
    
end
r_in = trans_lims(1:end-1);
r_ex = trans_lims(2:end);
ringvolume = binsize*pi*(r_ex.^2-r_in.^2);
chargematrix = fliplr(chargematrix);
densitymatrix = (chargematrix./ringvolume')*5e6;

tbin = linspace(trim_window_size,0,ceil(trim_window_size/binsize));
xlims = [tbin(1),tbin(end)];
ylims = [-trans_lims(end),trans_lims(end)];
P.plot_field_density('density_plot',densitymatrix,'r_plot',10*ylims,'z_plot',xlims);
P.plot_field_density('density_plot',chargematrix,'r_plot',10*ylims,'z_plot',xlims);

P.plot_name = '12dd3'; P.save_format = 'png'; P.fig_handle = gcf; 
% P.save_plot();
% save(['loading_files/',datadir,'densitytimeprofile.mat'],'densitymatrix','xlims','ylims')

charge = O.cylindrical_integration(trans_lims(1:end-1),fliplr(tbin),densitymatrix);

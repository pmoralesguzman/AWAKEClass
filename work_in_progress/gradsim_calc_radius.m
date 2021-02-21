%________________________________________________________________________
%
% Script to calculate the radius of the proton bunch.
%
% Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 18/02/2021
%________________________________________________________________________

% save directories
clear;
close all;

% paramaters
plasmaden = 1.81e14;
datadir = 'gp20';
dump_list = 0:3:123;
binsize = OsirisDenormalizer.ps2cm(0.41214);
trans_lims = linspace(0,0.16,74); %cm (74 for 0.16 cm (see above))
prop_distances = 9.958450708031729*dump_list;

r_in = trans_lims(1:end-1);
r_ex = trans_lims(2:end);
ringvolume = binsize*pi*(r_ex.^2 - r_in.^2);

xi_range = [18.5335,0];
dataformat = 'h5';


% get macroparticle info

O = OsirisDenormalizer('plasmaden',plasmaden,'datadir',datadir,...
    'property','raw','xi_range',xi_range,...
    'dump',0,'dataformat',dataformat,'species','proton_beam');
P = Plotty('plasmaden',plasmaden,'datadir',datadir,...
    'species','proton_beam',...
    'property_plot','density');


O.raw_dataset = 'x'; O.direction = 'r'; O.getdata(); O.assign_raw();
O.raw_dataset = 'x'; O.direction = 'z'; O.getdata(); O.assign_raw();
O.raw_dataset = 'p'; O.direction = 'r'; O.getdata(); O.assign_raw();
O.raw_dataset = 'q'; O.getdata(); O.assign_raw();
O.raw_dataset = 'ene'; O.getdata(); O.assign_raw();

% to get simulation window size, needed for the distance
O.dataformat = 'mat'; % done for gm20
O.property = 'density'; O.getdata(); O.trim_data(); O.denorm_distance();
O.dataformat = 'h5';


% manual trim of raw data
ixi = (O.denorm_distance(O.nz_raw) - O.dtime) < xi_range(1);
z = O.denorm_distance(O.nz_raw(ixi));
r = O.denorm_distance(O.nr_raw(ixi));
q = O.q_raw(ixi);

O.npr_raw = O.npr_raw(ixi);
O.nr_raw = O.nr_raw(ixi);
O.nE_raw = O.nE_raw(ixi);

    % delay_time = (O.dtime + O.simulation_window - z)/O.c_cm*1e12; % ps
    % t_simulation_window = (O.simulation_window)/O.c_cm*1e12; % ps
%     z_pos = (O.dtime + O.simulation_window - z); %z with respect to time
    z_pos = z - O.dtime;
    trim_window_size = max(O.z) - min(O.z);

for n = 1:length(dump_list)
        
    new_r = O.denorm_distance(O.charge_pusher(O,prop_distances(n)));
    rms_raw(n) = sqrt(sum((new_r.^2).*q)/sum(q))/sqrt(2);

%     chargematrix = zeros(length(trans_lims)-1,ceil(trim_window_size/binsize));
%     for rr = 1:(length(trans_lims)-1)
%         ir = (new_r >= trans_lims(rr)) & (new_r < trans_lims(rr+1));
%         if sum(ir) == 0; continue; end
%         
%         % select only those particles inside the transverse limits
%         q_r = q(ir);
%         % delay_time_r = delay_time(ir);
%         z_r = z_pos(ir);
%         
%         %[~,~,ind_bin] = histcounts(z_r,0:binsize:trim_window_size);
%         %A = accumarray(ind_bin+1,z_r);
%         % bin
%         ind_bin = ceil(z_r/binsize);
%         
%         % sort the indeces of the binning and then sort the
%         % elements in q
%         [ind_sort,ind_order] = sort(ind_bin);
%         q_sort = q_r(ind_order);
%         
%         % say how many there is in each bin in a cumulative way,
%         % to use that as indices
%         % for q, to quickly build up the indexed density along z
%         ind_sum = [0,0,histcounts(ind_sort,1:1:max(ind_bin),'Normalization','cumcount')];
%         
%         % go through each element in bz and arrange that charge
%         % according to the bins, using the indeces from histconts
%         for bz = 1:max(ind_bin)
%             chargematrix(rr,bz) = sum(q_sort(ind_sum(bz)+1:ind_sum(bz+1) ));
%         end
%         O.progress_dump('building density',rr,length(trans_lims)-1)
%         
%     end % trans lims (build profile)
% 
%     chargematrix = fliplr(chargematrix);
%     
%     densitymatrix = (chargematrix./ringvolume');
%     imagesc(chargematrix)
%     trans_profile = sum(densitymatrix,2);
%     trans_profile_fit = [flipud(trans_profile(:));trans_profile(:)];
%     x = [-fliplr(r_ex),r_ex];
%     f = fit(x(:),trans_profile_fit(:),'gauss1');
%     plot(x,trans_profile_fit);
%     rms_density(n) = f.c1/sqrt(2);
    
    
end


rms_theorie = 0.02*sqrt(1+(prop_distances/100).^2/4.9^2);



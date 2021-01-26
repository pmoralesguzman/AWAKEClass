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
measurement_point = 1350; % in cm, from plasma beginning
dump = 100;
binsize = 0.41214;
trans_lims = linspace(0,0.3,138); %cm (74 for 0.16 cm (see above))
xi_range = [18.5335,0];
dataformat = 'h5';

% get macroparticle info

O = OsirisDenormalizer('plasmaden',1.81e14,'datadir',datadir,...
    'property','raw','xi_range',xi_range,...
    'dump',dump,'dataformat',dataformat,'species','proton_beam');
P = Plotty('plasmaden',1.81e14,'plots_dir','slit','plot_name','3D','save_format','png');
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
O.nE_raw = O.nE_raw(ixi);

% front and back propagate

prop_distances = measurement_point + O.simulation_window - z;

new_r = O.denorm_distance(O.charge_pusher(O,prop_distances));

dy = (4000/2)*1e-4; % 74 um width

n_points = 10;

eq_spaced_v = linspace(0,1,n_points);
th_end_points = asin(dy./new_r);
i_complex = (imag(th_end_points) ~= 0);
th_end_points(i_complex) = pi/2;
q_temp = O.q_raw.*th_end_points;

th_mat = eq_spaced_v.*(th_end_points);
new_x = new_r.*cos(th_mat); 
% new_r = new_x(:);

new_y = new_r.*sin(th_mat); 

new_q = q_temp.*(ones(1,n_points)./n_points);
new_q = new_q(:);

new_z = z.*(ones(1,n_points));
% new_z = new_z(:);


delay_time = (O.dtime + O.simulation_window - new_z)/O.c_cm*1e12; % ps
t_simulation_window = (O.simulation_window)/O.c_cm*1e12; % ps

%% 
N = 15;
% random_picker = randperm(N,N); %length(O.q_raw)
random_picker = repmat(randperm(length(O.q_raw),N),n_points/2,1); random_picker = random_picker(:);
random_picker30 = repmat((1:2:n_points)',N,1); %randi(2,1,N);

p_x = new_x(sub2ind(size(new_x),random_picker,random_picker30));
p_y = new_y(sub2ind(size(new_y),random_picker,random_picker30));
p_z = new_z(sub2ind(size(new_z),random_picker,random_picker30));

cmp = colormap(parula(N));
cmp = repmat(cmp,n_points/2,1); 
close;

plot_x = [p_x;p_x;-p_x;-p_x];
plot_y = [p_y;-p_y;p_y;-p_y];
plot_z = [p_z;p_z;p_z;p_z];

plot_cmp = [cmp;cmp;cmp;cmp];



S = 50;
hold on
% plot3([min(plot_x),max(plot_x)],[0.0001 0],[0.0001 0],'LineWidth',3);

% scatter3sph(plot_z,plot_x,plot_y,'size',S,'color',plot_cmp)
scatter3(plot_z,plot_x,plot_y,S,plot_cmp,'filled')
% ylim([-0.05,0.05]);
% zlim([-0.05,0.05]);
% xlim();
hold off

P.fig_handle = gcf;
P.plot_name = 'Student3';
P.save_plot();


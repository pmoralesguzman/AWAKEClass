%________________________________________________________________________
% FFT of the proton or plasma electrons density distribution, or of a
% lineout of the wakefields.
%
% Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 08/06/2020
%________________________________________________________________________


% data directory
datadirs = {'gm20'};
plots_dir = ['charge_trans_evolution/',datadirs{1}];
P = Plotty();


% parameters
plasma_density = 1.81e14;

% properties
property = 'density';
species = 'proton_beam';

% simulation parameters
dump_list = 0:1:133;
dataformat = 'mat';
useAvg = false;

% limits
plasma_radius = 0.149; % cm
nslices = 201;
xi_range = [21 0.74];

% analysis parameters
scan_type = 'slice'; % slice, cumulative
on_axis = 'sum'; % int, sum, intw, lineout


% calculated parameters
trans_lims = 2*3.9598e-04:2*3.9598e-04:0.15;

% initialize variables
charge_in_slice = ones(length(trans_lims),length(dump_list));
prop_distance_m = zeros(1,length(dump_list));

AC = AwakeFFT(...
    'plasmaden',plasma_density,'property',property,'species',species,...
    'dataformat',dataformat,'useAvg',useAvg,...
    'trans_lims',trans_lims,'xi_range',xi_range,...
    'scan_type',scan_type,'on_axis',on_axis); 


for d = 1:length(datadirs)
    AC.datadir = datadirs{d};

    for n = 1:length(dump_list)
        AC.dump = dump_list(n);
        AC.fft_dataload();
        prop_distance_m(n) = AC.propagation_distance/100; % propagation distance in m

        switch AC.property
            case 'density'
                charge_in_slice(:,n) = AC.dz*sum(AC.fft_densitymatrix,2);
        end % switch property 
        AC.progress_dump('dump',n,length(dump_list));
    end %  for dump list
    
    fct = figure(1);
    [X,Y] = meshgrid(prop_distance_m,trans_lims*10);
    waterfall(X',Y',charge_in_slice');
    
    ylabel('trans. distance (mm)');
    xlabel('prop. distance (m)');
    zlabel('charge density (protons/cm)');
    
    drawnow;
    P.plots_dir = plots_dir;
    P.plot_name = ['charge_trans',AC.datadir,'n',num2str(dump_list(end))];
    P.fig_handle = fct;
    
%     P.save_plot();
    
end % for datadirs








%__________________________________________________________________________
% script to play around with the tracks data to create the laoding code
% (findings already implemented in the OsirisLoader class)
% For use with: Osiris 4.4.4
%
% P. I. Morales Guzman
% Last update: 22/06/2020
%
%__________________________________________________________________________


trackpath = 'gm20/MS/TRACKS/proton_beam-tracks-repacked1.h5';
% g0 particles: 47732
% save each ndumpfac/10 iterations
% total: 1000 points per particle

% g0: 1816365 908182
% gm20: 3070163 1535081
% gp20: 2647281 1323640

tracks_data_temp = h5read(trackpath,'/data');
tracks_iter_temp = h5read(trackpath,'/itermap');

%%
nsave = max(tracks_iter_temp(2,:)); % data was dumped each niter*nsave iterations (50)
itersave = unique(tracks_iter_temp(3,:)); % first iteration at which the dump was saved (1x20)
niter = itersave(2)/nsave;
ntracks_per_particle = length(itersave)*nsave; % track points per particle (1000)
npar = max(tracks_iter_temp(1,:));
iter_length = length(tracks_iter_temp);


tracks_in_order = zeros(npar,8,ntracks_per_particle);
ind = ones(npar,2);
begin = 1; final = 0;
ind(:,2) = zeros(npar,1);

indtable = ones(npar,2*ntracks_per_particle/nsave);
indtable(:,2) = zeros(npar,1);
par_occurance_table = zeros(npar,1);
% progressbar;

for ii = 1:iter_length
    
    par = tracks_iter_temp(1,ii);

    final = final + tracks_iter_temp(2,ii);
    ind(par,2) = ind(par,2) + tracks_iter_temp(2,ii);
    
    par_occurance_table(par) = par_occurance_table(par) + 1;
    
    indtable(par,par_occurance_table(par)*2+2) = ...
        indtable(par,(par_occurance_table(par)-1)*2+2) + tracks_iter_temp(2,ii);
    
   
    tracks_in_order(par,:,ind(par,1):ind(par,2))...
        = tracks_data_temp(:,begin:final);
    
    begin = final + 1;
    ind(par,1) = ind(par,2) + 1;
    
    
    indtable(par,par_occurance_table(par)*2+1) = ...
        indtable(par,(par_occurance_table(par)-1)*2+2) + 1;
    

    
%     progressbar(ii/iter_length)
    disp(ii)
end

tracks = tracks_in_order(:,5,:);
traplot = tracks;
traplot(traplot==0) = nan;
traplot = squeeze(traplot)';
plot(traplot(:,1:4))
% n = 1;
% while n <= length(trackiter)
%     if all(trackiter(2,))
% end




% time
% charge
% energy
% z
% r
% pz
% pr
% pth

%
% h = 0;
% for ii = 1:200
%    if sum(trackiter(1,:)==ii) == 20
%       h = [h,ii];
%    end
% end

% particle
%
% iteration number

% HDF5 proton_beam-tracks-repacked.h5
% Group '/'
%     Attributes:
%         'TYPE':  'tracks-2'
%         'NAME':  'proton_beam                                                     '
%         'NTRACKS':  47732
%         'NDUMP':  36928
%         'NITER':  3693
%         'DT':  0.007000
%         'XMIN':  0.000000 -0.005000
%         'XMAX':  532.000000 3.995000
%         'PERIODIC':  0 0
%         'MOVE C':  1 0
%         'QUANTS':  'n               ', 't               ', 'q               ', 'ene             ', 'x1              ', 'x2              ', 'p1              ', 'p2              ', 'p3              '
%         'UNITS':  '                ', '1/\omega_p      ', 'e               ', 'm_e c^2         ', 'c/\omega_p      ', 'c/\omega_p      ', 'm_e c           ', 'm_e c           ', 'm_e c           '
%     Dataset 'data'
%         Size:  8x38201057
%         MaxSize:  8x38201057
%         Datatype:   H5T_IEEE_F64LE (double)
%         ChunkSize:  []
%         Filters:  none
%         FillValue:  0.000000
%     Dataset 'itermap'
%         Size:  3x783470
%         MaxSize:  3x783470
%         Datatype:   H5T_STD_I32LE (int32)
%         ChunkSize:  []
%         Filters:  none
%         FillValue:  0
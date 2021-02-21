%________________________________________________________________________
% gradsim paper
% Script to follow the zero crossing of fields from waterfall plots. It
% also creates the waterfall plots if not yet made. 
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

O = OsirisDenormalizer('datadir','gm20','property','density',...
    'dump',0,'plasmaden',1.81e14);

O.getdata(); O.assign_density(); O.denorm_density(); O.denorm_distance();

O.property = 'raw'; O.dataformat = 'h5';
O.raw_dataset = 'x'; O.direction = 'r'; O.getdata(); O.assign_raw();
r = O.denorm_distance(O.nr_raw);
O.raw_dataset = 'q'; O.getdata(); O.assign_raw();
q = O.q_raw;
rms = sqrt(sum((r.^2).*q)/sum(q))/sqrt(2)


p = O.proton_beam;

pm = [flipud(p);p];

figure(1);
imagesc(pm);

% lo = pm(:,35300);
lo = sum(pm,2);
x = [-fliplr(O.r),O.r];

figure(2);
plot(x,lo)

f = fit(x(:),lo(:),'gauss1');

cr = O.radial_integration(O.r,O.z,p);
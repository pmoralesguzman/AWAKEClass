

O = OsirisDenormalizer('datadir','gm20','property','density',...
    'dump',0,'plasmaden',1.81e14);

O.getdata(); O.assign_density(); O.denorm_density(); O.denorm_distance();

p = O.proton_beam;

pm = [flipud(p);p];

figure(1);
imagesc(pm);

lo = pm(:,35300);
x = [-fliplr(O.r),O.r];

figure(2);
plot(x,lo)

f = fit(x(:),lo(:),'gauss1');

cr = O.radial_integration(O.r,O.z,p);
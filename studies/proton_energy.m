%%%%%%

datadir = 'gp20';
plasmaden = 1.81e14;
dump_list = 0:1:100;
bins = 21;
p_mass_GeV = 0.938;
% get initial protons

O = OsirisDenormalizer('plasmaden',plasmaden,'datadir',datadir,...
    'property','density',...
    'dump',0,'dataformat','mat','species','proton_beam');
O.getdata(); O.assign_density(); O.denorm_density(); O.denorm_distance();

P = Plotty('plasmaden',plasmaden,'plots_dir',['proton_energy/',datadir]);

initial_protons = O.cylindrical_integration(O.r,O.z,O.proton_beam);
original_edges = linspace(min(O.z),max(O.z),bins+1);
edges = original_edges;
edges(1) = -inf; edges(end) = inf;

O.property = 'raw'; O.raw_dataset = 'q'; O.dataformat = 'h5';
O.getdata(); O.assign_raw(); 

O.raw_dataset = 'ene';
O.getdata(); O.assign_raw();

initial_energy = sum(O.nE_raw.*O.q_raw)./sum(O.q_raw);
initial_protons_factor = initial_protons/sum(O.q_raw);

energy_per_proton = zeros(1,length(dump_list));

total_energy = zeros(1,length(dump_list));
z_pos = zeros(1,length(dump_list));

q_bin = zeros(length(dump_list),bins);
total_ene_bin = zeros(length(dump_list),bins);

switch datadir
    case 'R2gap0_2e14'
        title_string = ['Density Step (no gap)'];
    case 'R2gap100_2e14'
        title_string = ['Density Step (1 m gap)'];
    case 'R2gap0_2e14_nods'
        title_string = ['Constant density (no gap)'];
    case 'R2gap0g5'
        title_string = ['gradient = 0.5 %/m (no gap)'];
    otherwise
        title_string = '';
end

for n = 1:length(dump_list)
    
    O.dump = dump_list(n);
    O.raw_dataset = 'q';
    O.getdata(); O.assign_raw();
    
    O.raw_dataset = 'ene';
    O.getdata(); O.assign_raw();
    
    O.raw_dataset = 'x'; O.direction = 'z';
    O.getdata(); O.assign_raw();
    
    z = O.denorm_distance(O.nz_raw - O.ntime);
    O.denorm_distance();
    z_pos(n) = (O.dtime + O.simulation_window)/100; % m
    
    indz = ceil(z);
    unique_indz = unique(indz);
    
    
    energy_per_proton(n) = (sum((O.nE_raw + 1).*O.q_raw)./sum(O.q_raw))*p_mass_GeV; % GeV
    total_energy(n) = sum((O.nE_raw).*O.q_raw).*initial_protons_factor.*O.c_m^2*O.p_mass_kg;
    
    [~,~,ind_bin] = histcounts(z,edges);
    total_ene_bin(n,:) = accumarray(ind_bin,(O.nE_raw + 1).*O.q_raw);
    q_bin(n,:) = accumarray(ind_bin,O.q_raw);
    
    O.progress_dump('proton energy',n,length(dump_list));
end

ene_bin = (total_ene_bin./q_bin).*p_mass_GeV; % GeV

%% plot part 

fig_pE = figure;
plot(z_pos,energy_per_proton,'LineWidth',2);
xlim([min(z_pos),max(z_pos)])
xlabel('z');
ylabel('energy per proton (GeV/c^2)');
ylim([399.7,400.01]);
title(title_string);
P.plot_name = ['energy_proton'];
P.fig_handle = fig_pE;
P.save_plot();


fig_totene = figure;
plot(z_pos,total_energy,'LineWidth',2);
xlim([min(z_pos),max(z_pos)])
ylim([0,14000]);
xlabel('z');
ylabel('total kinetic energy (J)');
title(title_string);
P.plot_name = ['total_energy'];
P.fig_handle = fig_totene;
P.save_plot();


fig_watene = figure;
imagesc(original_edges(2:end),z_pos,fliplr(ene_bin));
set(gca,'YDir','normal');
set(gca,'XDir','reverse');
% xlim([min(original_edges),max(original_edges)])
colorbar_handle = colorbar;
ylabel('z (m)'); 
xlabel('\xi (cm)');
colorbar_handle.Label.String = ['energy per proton (GeV)'];
title(title_string);
P.plot_name = ['energy_waterfall'];
P.fig_handle = fig_watene;
P.save_plot();


fig_watch = figure;
imagesc(original_edges(2:end),z_pos,fliplr(q_bin.*initial_protons_factor));
caxis([0,max(q_bin.*initial_protons_factor,[],'all')])
colormap(flipud(gray));
set(gca,'YDir','normal');
set(gca,'XDir','reverse');
% xlim([min(original_edges),max(original_edges)])
colorbar_handle = colorbar;
ylabel('z (m)'); 
xlabel('\xi (cm)');
colorbar_handle.Label.String = ['protons'];
title(title_string);
P.plot_name = ['protons_waterfall'];
P.fig_handle = fig_watch;
P.save_plot();

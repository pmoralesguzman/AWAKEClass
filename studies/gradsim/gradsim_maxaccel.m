%---------------------------------------------------------------------
% Plot waterfall, energy gain from Ez along xi (point-like bunch),
% maximum energy gain and position, and phase evolution.
%
% Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 24/03/2020
%---------------------------------------------------------------------
clear;

load('loading_files/physical_constants.mat')

cc = [118,42,131
153,112,171
194,165,207
231,212,232
10,10,10
217,240,211
166,219,160
90,174,97
27,120,55]/256;

%% input parameters
dir_names = {'gm20','gm15','gm10','gm5',...
        'g0','gp5','gp10','gp15','gp20'};
    
% dir_names = {'gp5','gp10'};

plasma_density = 1.81e14;
dump_list = 0:100;
use_avg = false;
data_format = 'mat';
leg = {'-2 %/m','-1.5 %/m','-1 %/m','-0.5 %/m','0 %/m','0.5 %/m','1 %/m','1.5 %/m','2 %/m'};

%% calculated parameters
plasma_frequency = n2freq(plasma_density);
plasma_wavelength = 2*pi*c_cm/plasma_frequency;


for d = 1:length(dir_names)
dir_name = dir_names{d};

switch dir_name
    case  {'gm20'}
        title_gap = '(g = -20 %)';
    case 'gm15'
        title_gap = '(g = -15 %)';
    case 'gm10'
        title_gap = '(g = -9 %)';
    case 'gm5'
        title_gap = '(g = -5 %)';
    case  'g0'
        title_gap = '(g = 0 %)';
    case 'gp5'
        title_gap = '(g = +5 %';
    case 'gp9'
        title_gap = '(g = +9 %)';
    case 'gp15'
        title_gap = '(g = +13 %)';
    case 'gp20'
        title_gap = '(g = +20 %)';
end

%% waterfall plot
[waterfall_build,x_axis,L,im] = plot_waterfall('dir_name',dir_name,...
    'plasma_density',plasma_density,...
    'dump_list',dump_list,'use_avg',use_avg,...
    'data_format',data_format);
% title(['Waterfall plot ',title_gap])
% savefig(['plots/energy_gain/waterfall',dir_name,'.fig'])

x_copy = x_axis;
indxi = (x_axis > 13) & (x_axis < 15);
waterfall_build = waterfall_build(:,indxi);
x_axis = x_axis(indxi);

%% field integral plot
field_integral = cumtrapz(L/100,-(waterfall_build));
% figure;
% imagesc(field_integral)
% colormap(bluewhitered);
% ylabel('MV');
% xlabel('distance from plasma end (m)');
% title('Energy gain from the end to the beginning of the plasma')

%% max energy gain plot
[max_energygain,pos_maxenergygain] = max(field_integral);
% 
% figure(999);
% colororder(cc)
% hold all
% plot(x_axis,max_energygain,'LineWidth',2);
% hold off
% legend(leg,'location','best');
% title('Maximum energy gain along the bunch')
% ylabel('MeV')
% xlabel(' (cm)')
% savefig(['plots/energy_gain/maxene_alongbunch','.fig',])
% saveas(gcf,['plots/energy_gain/maxene_alongbunch','.png']);

[maxmax_energygain,best_position] = max(max_energygain);

% figure;
% plot(x_axis,L(pos_maxenergygain),'LineWidth',2)
% title('Position for maximum energy gain along the plasma')

% plot energy gain
figure(1000);
colororder(cc)
hold on
field_integral_save(d,:) = field_integral(:,best_position);
plot(L,field_integral_save(d,:),'LineWidth',2)
hold off
title('Energy gain for each injection position along the plasma')
legend(leg,'location','best');
ylabel('MeV');
xlabel('distance from plasma end (m)');

savefig(['collabmeeting/energy_gain/maxene_alongplasma','.fig',])
saveas(gcf,['collabmeeting/energy_gain/maxene_alongplasma','.png']);


%%_________ energy gain
% Also do for all xi, and create waterfall plot




%% field amplitude
figure(1001);
colororder(cc)
hold on
waterfall_build_save(d,:) = -flip(waterfall_build(:,best_position));
plot(L,waterfall_build_save(d,:),'LineWidth',2)
hold off
title('Ez along L for the max. energy gain')
legend(leg,'location','best');
xlabel('propagation length (m)')
ylabel('Ez MV/m')

savefig(['collabmeeting/energy_gain/field_alongplasma','.fig',])
saveas(gcf,['collabmeeting/energy_gain/field_alongplasma','.png']);


%% Find max amplitude around the best position

amplitude_indices = (x_axis > (x_axis(best_position) - plasma_wavelength/2)) & (x_axis < (x_axis(best_position) + plasma_wavelength/2));
waterfall_around_max = waterfall_build(:,amplitude_indices);
[max_amplitude,maxamp_ind] = max(-waterfall_around_max,[],2);

figure(1002)
colororder(cc)

hold on
max_amplitude_save(d,:) = flipud(max_amplitude);
plot(L,max_amplitude_save(d,:),'LineWidth',2);
hold off
legend(leg,'location','best');
title('Max. amplitude along the plasma close to injection point')
xlabel('propagation length (m)')
ylabel('Ez MV/m')

savefig(['collabmeeting/energy_gain/field_alongplasma','.fig',])
saveas(gcf,['collabmeeting/energy_gain/field_alongplasma','.png']);



%% __________ phase
% get min(abs(of the waterfall plot))
% diff_maxamplitude = diff(x_axis(max_amplitude_position));
% diff_maxamplitude_ind = find(diff_maxamplitude > 2*plasma_wavelength);
%


if true


maxamp_pos = fliplr(x_axis(maxamp_ind)-x_axis(1));
maxamp_pos = maxamp_pos(5:end);
diff_maxamp = diff(maxamp_pos);
diff_maxampplus_ind = find((diff_maxamp) > 0.8*plasma_wavelength);
diff_maxampminus_ind = find((diff_maxamp) < -0.8*plasma_wavelength);

for n = 1:length(diff_maxampplus_ind)
    maxamp_pos(1:diff_maxampplus_ind(n)) = maxamp_pos(1:diff_maxampplus_ind(n)) + plasma_wavelength;
end

for n = 1:length(diff_maxampminus_ind)
    maxamp_pos(1:diff_maxampminus_ind(n)) = maxamp_pos(1:diff_maxampminus_ind(n)) - plasma_wavelength;
end

L_maxamp = L(5:end);

figure(1005);
colororder(cc)

N = 12;

hold on
plot(L_maxamp,maxamp_pos,'LineWidth',2);
hold off
legend(leg,'location','best');
title('Phase evolution at max. amplitude')
xlabel('propagation length (m)')
ylabel('phase (rad)')
end

end








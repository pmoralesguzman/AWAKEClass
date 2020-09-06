%__________________________________________________________________________
% Example for the classes in the AWAKE Osiris Analysis Matlab Package
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

% 
datadir = 'benchmark_etc';
plasmaden = 7e14;

dump_list = 1:1:2;
useAvg = false;
dataformat = 'h5';

% simulation parameters
species = 'electrons';

% tag conditions
trans_range = [0 10];
% 0.01 c/wp at 7e14 = 2.0085e-04 cm
% for benchmark_etc, simulation window = 597.5 c/wp = 12.0010
% particles taken from 595.81 to 596.01 = 11.9671 to 11.9711 cm
xi_range = [11.9711-2.0085e-4/2.2,12]; %this is z - ct, not ct - z (meaning the numbers start at the back of the bunch)
ntag = 10000;

OPT = OsirisParticleTracking('datadir',datadir,'plasmaden',plasmaden,...
    'species',species,...
    'dataformat',dataformat,...
    'dump_list',dump_list,'useAvg',useAvg,...
    'trans_range',trans_range,'xi_range',xi_range,...
    'property','raw','ntag',ntag);

tags = zeros(OPT.ntag*(length(dump_list)),2);
temp_tags = cell(1,length(dump_list));

% loop for transverse segments in one dump
% for n = 1:length(dump_list)
%     for tr = 1:15
%         OPT.trans_range = trans_range;
%         OPT.dump = dump_list(n);
%         OPT.select_tags();
%         temp_tags{n,tr} = OPT.selected_tags';
%         tags((tr-1)*OPT.ntag+[1:OPT.ntag],:) = temp_tags{n,tr};
%         trans_range = [trans_range(2) trans_range(2)+0.01];
%         disp(tr)
%     end
% end

% loop for all particles in each dump
initial = 1; sum_tags = 0;
for n = 1:length(dump_list)
        OPT.dump = dump_list(n);
        OPT.ntag = ntag;
        OPT.select_tags();
        temp_tags{n} = OPT.selected_tags';
        final = initial + OPT.ntag - 1;
        tags(initial:final,:) = temp_tags{n};
        initial = final + 1;
        sum_tags = sum_tags + length(OPT.selected_tags);
        OPT.progress_dump('selecting tags',n,length(dump_list));
end

tags(tags(:,1) == 0, :) = [];

number_of_tags = length(tags);

% writematrix(['save_files/',datadir,'/electrons',datadir,'.tags'],[number_of_tags,0;tags],'delimeter', ',' ,'precision','%i')
writematrix([number_of_tags,0;tags],['save_files/',datadir,'/electrons',datadir,'t1m.tags'],'FileType','text')

%% check by plotting

% loop for transverse segments in one dump
% 
% r_total = 0;
% for n = 1:length(dump_list)
%     for tr = 1:15
%     OPT.in_taglist = temp_tags{n,tr};
%     OPT.dump = dump_list(n);
%     OPT.raw_dataset = 'x'; OPT.direction = 'r';
%     OPT.getdata(); OPT.assign_raw();
%     OPT.find_tags();
%     disp(n)
%     r = OPT.denorm_distance(OPT.nr_raw(OPT.ind_tag));
%     r_total = [r;r_total];
%     z = dump_list(n)*ones(length(r),1);
%     
%     hold on
%     scatter(z,r)
%     hold off
%     end
% end
% figure
% histogram(r_total);


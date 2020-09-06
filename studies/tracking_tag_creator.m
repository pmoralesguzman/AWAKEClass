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


datadir = 'gp20';
plasmaden = 1.81e14;
dump_list = 50:100;
reference_dump = 100;
useAvg = false;
dataformat = 'mat';
trans_range = [0.14 0.15];
ntag = 1000;

OPT = OsirisParticleTracking('datadir',datadir,'plasmaden',plasmaden,...
    'dataformat',dataformat,...
    'dump_list',dump_list,'reference_dump',reference_dump,'useAvg',useAvg,...
    'trans_range',trans_range,...
    'property','raw','ntag',ntag);

tags = zeros(ntag*(length(dump_list)+15),2);
for n = 1:length(dump_list)
    OPT.reference_dump = dump_list(n);
    if dump_list(n) == reference_dump(end)
        OPT.trans_range = [0 0.15];
    end
    OPT.select_tags();
    temp_tags{n} = OPT.selected_tags';
    tags((n-1)*ntag+[1:ntag],:) = temp_tags{n};
    disp(n)
end

tags(tags(:,1) == 0,:) = [];

tags = unique(tags,'rows');

number_of_tags = length(tags);

csvwrite('protonsgp20.tags',[number_of_tags,0;tags])

%% check by plotting


for n = 1:length(dump_list)
    OPT.in_taglist = temp_tags{n};
    OPT.dump = dump_list(n);
    OPT.raw_dataset = 'x'; OPT.direction = 'r';
    OPT.getdata(); OPT.assign_raw();
    OPT.find_tags();
    disp(n)
    r = OPT.denorm_distance(OPT.nr_raw(OPT.ind_tag));
    z = dump_list(n)*ones(length(r),1);
    
    hold on 
    scatter(z,r)
    hold off  
end


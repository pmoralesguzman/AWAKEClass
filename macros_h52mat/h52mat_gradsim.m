

%%%%%%
% Comments comments
%%%%



% datadirs = {'gm20_m','g0zh','g0r2','g0z2','g0dt92'}; % (27/11/2020) and all raw electrons from all convergence except g0rh 
% missing RAW for gm20_m
% datadirs = {'gm15','gm10','gm5','g0','gp5','gp10','gp15','gp20'};
datadirs = {'gm20'};
% properties = {'fields','density','raw'};
properties = {'raw'};
speciess = {'proton_beam'};
dump_list = 0:1:134;


O = OsirisLoader('plasmaden',1.81e14,'dataformat','h5');


for ii = 1:length(datadirs)
    % specify the name of the directories where the data is in cell form
    % ex.: {'AW_ionizationfront18HD','AW_onlyproton','AW_apseed_HD256','AW_apseedx5_HD256','AW_apseed20_HD256'}
    
    O.datadir = datadirs{ii};
    
    for jj = 1:length(properties)
        
        O.property = properties{jj};
        switch O.property
            case 'density'
                for kk = 1:length(speciess)
                    O.species = speciess{kk};
                    
                    for n = 1:length(dump_list)
                        O.dump = dump_list(n);
                        O.h52mat();
                        O.progress_dump('h52mat density',n,length(dump_list))
                    end
                end
            case 'fields'
                for kk = 1:3
                    switch kk
                        case 1; O.field = 'e'; O.direction = 'z';
                        case 2; O.field = 'e'; O.direction = 'r';
                        case 3; O.field = 'b'; O.direction = '\theta';
                    end
                    
                    for n = 1:length(dump_list)
                        O.dump = dump_list(n);
                        O.h52mat();
                        O.progress_dump('h52mat fields',n,length(dump_list))
                    end % n dumps
                end % kk fields
            case 'raw'
                for ll = 1:length(speciess)
                    O.species = speciess{ll};
                    for kk = 1:7
                        switch kk
                            case 1; O.raw_dataset = 'q';
                            case 2; O.raw_dataset = 'x'; O.direction = 'z';
                            case 3; O.raw_dataset = 'p'; O.direction = 'z';
                            case 8; O.raw_dataset = 'tag';
                            case 5; O.raw_dataset = 'ene';
                            case 6; O.raw_dataset = 'x'; O.direction = 'r';
                            case 7; O.raw_dataset = 'p'; O.direction = 'r';
                            case 4; O.raw_dataset = 'p'; O.direction = '\theta';
                        end
                        
                        for n = 1:length(dump_list)
                            O.dump = dump_list(n);
                            O.h52mat();
                            O.progress_dump('h52mat raw',n,length(dump_list))
                        end % n dumps
                    end %% ll speciess
                end % kk raw
        end
    end % length properties
    
end % length datadirs


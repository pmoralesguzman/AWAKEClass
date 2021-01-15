%%%%%%


O = OsirisLoader('datadir','R2gap0_2e14',...
    'property','raw','raw_dataset','q',...
    'dump',110,'dataformat','h5','species','electron_bunch');
O.getdata(); O.assign_raw();

O.raw_dataset = 'ene';

O.getdata(); O.assign_raw();

initial_energy = sum(O.nE_raw.*O.q_raw)./sum(O.q_raw);

x = 1;
for n = 110:200
    O.dump = n;
    O.raw_dataset = 'q';
    O.getdata(); O.assign_raw();
    
    O.raw_dataset = 'ene';
    O.getdata(); O.assign_raw();
    
%     energy(x) = sum(O.nE_raw.*O.q_raw)./sum(O.q_raw)/initial_energy;
    energy(x) = max(O.nE_raw)/initial_energy;
    x = x + 1;
    O.progress_dump('hola',n,200);
end

    plot(linspace(11,20,length(energy)),energy,'LineWidth',2);

    xlabel('z');
    ylabel('max energy / average initial energy (a.u.)');
    title('gap = 0 m')
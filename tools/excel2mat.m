%__________________________________________________________________________
% Little script to transform excel to mat files
% For use with: Osiris 4.4.4
%
% AWAKE Experiment
%
% Work in progress
%
% P. I. Morales Guzman
% Last update: 04/06/2020
%__________________________________________________________________________


excel_filename = 'rscan_gradients.xlsx';
excel_range = 'A3:K83';

table_contents = readtable(excel_filename,'range',excel_range,...
    'ReadVariableNames',true,'PreserveVariableNames',true);

table_contents(1:2,:) = []; 

radius_mm = table_contents.Var2;
radius_sigma = table_contents.Radius;

gp20 = table_contents.("0.02");
gp13 = table_contents.("0.013");
gp9 = table_contents.("0.009");
gp4 = table_contents.("0.004");
g0 = table_contents.("0");
gm5 = table_contents.("-0.005");
gm9 = table_contents.("-0.009");
gm15 = zeros(length(table_contents.("0")),1);
gm19 = table_contents.("-0.019");

save('loading_files/fvsgr_exp.mat','radius_mm','radius_sigma','gp20','gp13',...
    'gp9','gp4','g0','gm5','gm9','gm15','gm19');


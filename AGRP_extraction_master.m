%AGRP extraction master V2
% Code written by Christian Pedersen
% Code adapted by Jennifer Dezet Deem and Sarah Larsen

clear all
close all
clc
fclose('all')


%% AGRP thermo session 1, 22 30 10 repeat ramp, different order each day but same ramp

blockname1{1} = 'Agrp6-5-191115-104205';

filepath1 = 'C:\Users\Katie\Documents\MATLAB\Jennifer AGRP May 2020 GitHub';

%%

% specify order of polynomial for baseline drift fit (0 = horizontal line)
polyOrder = 2;


for nn = [1]

    savename1 = ['raw_photom_',num2str(nn),'.mat'];

    tdtExtract_fxn_v1(filepath1,blockname1{nn},savename1,polyOrder)

end















date = '20171207';

folder = '/home/alvita/Dropbox (Vision Lab Cal Tech)/yuli-pharma/';
% folder = 'D:/Dropbox/Dropbox (Vision Lab Cal Tech)/yuli-pharma/';
addpath([folder 'mcd/']);
folder_results = [folder 'results/' date '/'];

switch date
    case '20171207'
        load([folder 'MEA_data/20171207_ethanol/allspk_kilosorted_20171207.mat']);
        trials_cntrl = 5:9;
        trials_drug = 10:22;
        trials_washout = 29:38;
    case '20181215'
        load([folder 'MEA_data/20181215_ethanol/20181215_allspkbycell.mat']);
end

numcells = size(allspk_bycell, 1);
T = 38*29/60;
binsize = 0.1;
timepoints = ceil(T/binsize);

PSTH = {};
PSTH.binsize = binsize;
PSTH.cntrl = zeros(numcells, timepoints);
PSTH.drug = zeros(numcells, timepoints);
PSTH.washout = zeros(numcells, timepoints);

for cellnum = 1:numcells
    disp(cellnum);
    
    [psth_cntrl, edges, std_cntrl] = calc_psth(folder, allspk_bycell, cellnum, trials_cntrl, binsize, T);
    [psth_drug, edges, std_drug] = calc_psth(folder, allspk_bycell, cellnum, trials_drug, binsize, T);
    [psth_washout, edges, std_washout] = calc_psth(folder, allspk_bycell, cellnum, trials_washout, binsize, T);
    
    PSTH.cntrl(cellnum, :) = psth_cntrl;
    PSTH.drug(cellnum, :) = psth_drug;
    PSTH.washout(cellnum, :) = psth_washout;
    
    PSTH.cntrl_std(cellnum, :) = std_cntrl;
    PSTH.drug_std(cellnum, :) = std_drug;
    PSTH.washout_std(cellnum, :) = std_washout;
end
PSTH.edges = edges;

save([folder_results 'PSTH.mat'], 'PSTH');
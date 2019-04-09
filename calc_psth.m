function [psth, edges, variability] = calc_psth(folder, allspk_bycell, cellnum, trials, binsize, T)

if nargin == 4
    binsize = [];
    psth = [];
else
    edges = 0:binsize:T+binsize; 
    psth = zeros(length(trials), length(edges)-1);
end

for t=1:length(trials) 
    trial = trials(t);

    %% get sync data for trial
    syncpath = [folder 'MEA_data/20171207_ethanol/20171207_processed/data' num2str(trial)];
    sync = LoadConvertedData(syncpath);
    sync = GetStimTiming_yulimod(sync);

    %% align spikes to stimulus
    spikes = allspk_bycell(cellnum, trial);
    spikes = spikes{1} - sync.triggers(1)/10000;

    %% bin spikes
    if isempty(binsize)
        edges = sync.triggers/10000 - sync.triggers(1)/10000;
        binsize = mean(diff(sync.triggers)) / 10000;
    end
    psth(t, :) = histcounts(spikes, edges)/binsize;
    
end

variability = std(psth, 0, 1);
psth = mean(psth, 1);

end
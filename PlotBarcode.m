folder = '/home/alvita/Dropbox (Vision Lab Cal Tech)/yuli-pharma/';
load([folder 'MEA_data/20171207_ethanol/Param_kilosort_20171207.mat'])
load([folder 'MEA_data/20171207_ethanol/allspk_kilosorted_20171207.mat'])

row = 2;
for cellnum = 1
    figure;
    
    for trial=1:length(Param.DataInfo)

%         %% load parameter file for trial
%         load([folder 'param_stimuli/20171207/Checkerboard_' str(trial) '.mat'])
%         stimparam = param;

        if ~isempty(strfind(Param.DataInfo{1, trial}, 'barcode'))
            
            %% get sync data for trial
            syncpath = [folder 'MEA_data/20171207_ethanol/20171207_processed/data' num2str(trial)];
            sync = LoadConvertedData(syncpath);
            sync = GetStimTiming_yulimod(sync);
            
            %% align spikes to stimulus
            spikes = allspk_bycell(cellnum, trial);
            spikes = spikes{1} - min(sync.triggers(1)/10000);
            
            %% plot spikes
            plot_raster(spikes, row)
            row = row+1;
        end
    end
    
    for row=[6, 22, 34, 47]
        plot([0 20], [row row], 'm', 'LineWidth', 3)
    end
    xlim([0 20])

end
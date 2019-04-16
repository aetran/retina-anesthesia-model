folder = '/home/alvita/Dropbox (Vision Lab Cal Tech)/yuli-pharma/';
% folder = 'D:/Dropbox/Dropbox (Vision Lab Cal Tech)/yuli-pharma/';

%% make model choices
% method_denoise: parametric, svd
% method_nonlinear: sigmoid, relu
% method_stimfitnl: flicker, barcode
date = '20171207';
method_denoise = 'parametric';
method_nonlinear = 'sigmoid';
method_stimfitnl = 'barcode';
if strcmp(method_denoise, 'svd')
    SVD = struct('ntimepoints_to_keep', 16, 'npixels_to_keep', 10, 'components_to_use', [1, 2]);
else
    SVD = [];
end
opt_figure = 0;
opt_close = 1;
opt_save = 1;
cells = [1:27 29:153];

%% load STAs and reconstructed flicker bar stimulus and barcode stimulus
addpath([folder 'mcd/']);
fps = 60;
T = 38*29/60; % 20 seconds at 60 fps
xq = 1/fps:1/fps:T;
    
switch date
    case '20171207'
        load([folder 'MEA_data/20171207_ethanol/allspk_kilosorted_20171207.mat']);
        trial_sta = 4;
    case '20181215'
        load([folder 'MEA_data/20181215_ethanol/20181215_allspkbycell.mat']);
end
folder_results = [folder 'results/' date '/'];
load([folder_results 'ST.mat']);
load([folder_results 'stim.mat']);
load([folder_results 'PSTH.mat']);
folder_figresults = [folder_results 'figs-LN_' method_denoise '-' method_nonlinear '-' method_stimfitnl '/'];
if ~exist(folder_figresults, 'dir')
    mkdir(folder_figresults);
end

%% figure settings
set(groot,'defaultLineLineWidth',2)
col = ones(65,3); 
col(1:33,2)=0:1/32:1; 
col(1:33,3)=0:1/32:1;
col(33:65,1) = 1:-1/32:0; 
col(33:65,2) = 1:-1/32:0;
col = flipud(col);

%% loop over cells
numcells = ST.Ncell;
Nx = ST.Nx;
Nw = ST.Nw;
goodnessoffit = zeros(numcells, 1);
paramsspace = zeros(Nx, numcells);
paramstime = zeros(numcells, ST.Nw);
paramsl = zeros(numcells, ST.Nx, ST.Nw);
switch method_nonlinear
    case 'relu'
        paramsnl = zeros(numcells, 2);
    case 'sigmoid'
        paramsnl = zeros(numcells, 3);
end
if isempty(cells)
    cells = 1:numcells;
end

for cellnum = cells
    disp(cellnum);
    
    psth_cntrl_interp = interpolate_signal(PSTH.cntrl(cellnum, :), PSTH.edges, xq);
    psth_drug_interp = interpolate_signal(PSTH.drug(cellnum, :), PSTH.edges, xq);
    
    %% L: get STA of specific cell
    switch method_denoise
        case 'parametric'
            [raw_space, raw_time, profile_space, profile_time, profile_spacetime] = sta_parametric(cellnum, ST);
            [~, output_linear_fitnl] = linear_filter(profile_space, profile_time, stimulus_flckbr);
            [~, output_linear] = linear_filter(profile_space, profile_time, stimulus_barcode);
            paramsspace(:, cellnum) = profile_space;
            paramstime(cellnum, :) = profile_time;
        case 'svd'
            profile_spacetime = sta_svd(cellnum, ST, SVD.ntimepoints_to_keep, SVD.npixels_to_keep, SVD.components_to_use);
            output_linear_fitnl = conv_sta(profile_spacetime, stimulus_flckbr);
            output_linear = conv_sta(profile_spacetime, stimulus_barcode);
    end
    paramsl(cellnum, :, :) = profile_spacetime;

    %% N: fit nonlinearity on flicker bar
    switch method_stimfitnl
        case 'barcode'
            output_linear_fitnl = output_linear;
%             psth_fitnl = psth_cntrl_interp;
            trials_cntrl = 5:9;
            binsize = 1/fps;
            [psth_fitnl, ~, ~] = calc_psth(folder, allspk_bycell, cellnum, trials_cntrl, binsize, T);
        case 'flicker'
            [psth_fitnl, ~, ~] = calc_psth(folder, allspk_bycell, cellnum, trial_sta);
    end
    [fun_nonlinear, pi_nonlinear] = nonlinearity_fit(method_nonlinear, output_linear_fitnl, psth_fitnl);
    paramsnl(cellnum, :) = pi_nonlinear;

    %% predicted and actual firing rates
    output_nonlinear = fun_nonlinear(pi_nonlinear, output_linear);
   
    rho = corrcoef(output_nonlinear, psth_cntrl_interp);
    goodnessoffit(cellnum) = rho(1, 2);

    %% figure
    if opt_figure
        figure('units','normalized','outerposition',[0 0 1 1]);
        
        % firing rates
        subplot(3, 1, 2:3); 
        hold on;
        legendtext = sprintf('LN: R = %0.2f', rho(1, 2));
        plot(xq, output_nonlinear, 'k', 'DisplayName', legendtext);
        plot(xq, psth_cntrl_interp, 'b', 'DisplayName', 'control');
        plot(xq, psth_drug_interp, 'r', 'DisplayName', 'drug');
        xlim([0, T-0.1]);
        legend('show');
        xlabel('Time (s)');
        ylabel('Firing Rate (Hz)');

        % STA
        subplot(3, 5, 1);
        sta_vec = ST.average{2};
        spatial_temporal_map = reshape(sta_vec(:, cellnum)/ max(sta_vec(:, cellnum)), ST.Nx, ST.Nw);
        imagesc(spatial_temporal_map(:, :));
        colormap(col);
        caxis([-1,1]);
        xlabel('Frames')
        ylabel('Bar Position');
        yticklabels({});
        
        if strcmp(method_denoise, 'parametric')
            subplot(3, 5, 2); hold on;
            plot(raw_space, 1:numel(profile_space), 'b');
            plot(profile_space, 1:numel(profile_space), 'k');
            ylabel('Bar Position');
            yticklabels({});
            ylim([1, numel(profile_space)]);
            set(gca, 'YDir', 'reverse')

            subplot(3, 5, 3); hold on;
            plot(0:numel(profile_time)-1, raw_time, 'b');
            plot(0:numel(profile_time)-1, profile_time, 'k');
            xlabel('Frames');
            xlim([0, numel(profile_time)-1]);
        end
        
        subplot(3, 5, 4);
        imagesc(profile_spacetime);
        colormap(col);
        caxis([-1, 1]);
        xlabel('Frames')
        ylabel('Bar Position');
        yticklabels({});
        
        % nonlinearity
        subplot(3, 5, 5); hold on;
        scatter(output_linear_fitnl, psth_fitnl, 'MarkerFaceColor','b','MarkerEdgeColor','b', 'MarkerFaceAlpha',.2, 'MarkerEdgeAlpha',.2);
        input_range = min(output_linear_fitnl):0.1:max(output_linear_fitnl);
        output_nonlinear_flckbr = fun_nonlinear(pi_nonlinear, input_range);
        plot(input_range, output_nonlinear_flckbr, 'k');
        xlabel('Linear Prediction');
        ylabel('Firing Rate (Hz)');

        % finish and save
        suptitle(num2str(cellnum));
        if opt_close
            saveas(gcf, [folder_figresults, num2str(cellnum), '.jpg']);
            close;
        end
    end
end
if opt_save
    save([folder_results 'LN_' method_denoise '-' method_nonlinear '-' method_stimfitnl '.mat'], 'paramsnl', 'paramsl', 'paramsspace', 'paramstime', 'goodnessoffit', 'SVD');
end
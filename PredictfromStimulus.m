%% make model choices
% method_denoise: parametric, svd
% method_nonlinear: sigmoid, relu

method_denoise = 'parametric';
method_nonlinear = 'sigmoid';
SVD = struct('ntimepoints_to_keep', 16, 'npixels_to_keep', 10, 'components_to_use', [1, 2]);

opt_figure = 0;
opt_close = 1;

%% get spike data from barcode stimulus
folder = '/home/alvita/Dropbox (Vision Lab Cal Tech)/yuli-pharma/';
% folder = 'D:/Dropbox/Dropbox (Vision Lab Cal Tech)/yuli-pharma/';
load([folder 'MEA_data/20171207_ethanol/allspk_kilosorted_20171207.mat']);
fps = 60;
T = 38*29/60*fps; % 20 seconds at 60 fps

%% load reconstructed flicker bar stimulus and barcode stimulus
load('20171207_ethanol_stimulus.mat')

%% get ST from flicker bar stimulus
load('20171207_ethanol_ST.mat')

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
cells = 1:numcells;

for cellnum = cells
    disp(cellnum);
    
    %% L: get STA of specific cell
    switch method_denoise
        case 'parametric'
            [raw_space, raw_time, profile_space, profile_time, profile_spacetime] = sta_parametric(cellnum, ST);
            [~, output_linear_flckbr] = linear_filter(profile_space, profile_time, stimulus_flckbr);
            [~, output_linear] = linear_filter(profile_space, profile_time, stimulus_barcode);
            paramsspace(:, cellnum) = profile_space;
            paramstime(cellnum, :) = profile_time;
        case 'svd'
            profile_spacetime = sta_svd(cellnum, ST, SVD.ntimepoints_to_keep, SVD.npixels_to_keep, SVD.components_to_use);
            output_linear_flckbr = conv_sta(profile_spacetime, stimulus_flckbr);
            output_linear = conv_sta(profile_spacetime, stimulus_barcode);
    end
    paramsl(cellnum, :, :) = profile_spacetime;

    %% N: fit nonlinearity on flicker bar
    trial = 4;
    [psth_flckbr, ~, ~] = PSTH_aet(folder, allspk_bycell, cellnum, trial);
    [fun_nonlinear, pi_nonlinear] = FitNonlinearity_aet(method_nonlinear, output_linear_flckbr, psth_flckbr);
    paramsnl(cellnum, :) = pi_nonlinear;

    %% LN output: control
    output_nonlinear = fun_nonlinear(pi_nonlinear, output_linear);
    xq = 1/fps * (1:T);

    %% control trials
    binsize = 0.1;
    trials = 5:9;
    [psth_cntrl, edges, rho_cntrl, std_cntrl] = PSTH_aet(folder, allspk_bycell, cellnum, trials, binsize, T);    
    psth_cntrl_interp = interpolate_signal(psth_cntrl, edges, xq);
    std_cntrl_interp = interpolate_signal(std_cntrl, edges, xq);

    %% drug trials
    trials = 10:22;
    [psth_drug, edges, rho_drug, std_drug] = PSTH_aet(folder, allspk_bycell, cellnum, trials, binsize, T);
    psth_drug_interp = interpolate_signal(psth_drug, edges, xq);
    
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
        xlim([0, T/fps]);
        legend('show');
        xlabel('Time (s)');
        ylabel('Firing Rate (Hz)');

        % STA
        subplot(3, 4, 1);
        sta_vec = ST.average{2}; % kyu version only calc one channel    
        spatial_temporal_map = reshape(sta_vec(:, cellnum)/ max(sta_vec(:, cellnum)), ST.Nx, ST.Nw);
        imagesc(spatial_temporal_map(:, :));
        colormap(col);
        caxis([-1,1]); % set color bar the fixed range of -1 to 1 across trials
        xlabel('Frames')
        ylabel('Bar Position');
        yticklabels({});
        
        if strcmp(method_denoise, 'parametric')
            subplot(3, 4, 2); hold on;
            plot(raw_space, 1:numel(profile_space), 'b');
            plot(profile_space, 1:numel(profile_space), 'k');
            ylabel('Bar Position');
            yticklabels({});
            ylim([1, numel(profile_space)]);
            set(gca, 'YDir', 'reverse')

            subplot(3, 4, 3); hold on;
            plot(0:numel(profile_time)-1, raw_time, 'b');
            plot(0:numel(profile_time)-1, profile_time, 'k');
            xlabel('Frames');
            xlim([0, numel(profile_time)-1]);
        end
        
        subplot(3, 4, 4);
        imagesc(profile_spacetime);
        colormap(col);
        caxis([-1, 1]);
        xlabel('Frames')
        %     ylabel('Bar position 100um/bar')
        ylabel('Bar Position');
        yticklabels({});
        
        suptitle(num2str(cellnum));
        
        % finish and save
        if opt_close
            saveas(gcf, [num2str(cellnum), '.jpg']);
            close;
        end
    end
end
save(['20171207_LNfitparams_' method_denoise '_' method_nonlinear '.mat'], 'paramsnl', 'paramsl', 'paramsspace', 'paramstime', 'goodnessoffit');
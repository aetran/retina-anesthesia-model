folder = '/home/alvita/Dropbox (Vision Lab Cal Tech)/yuli-pharma/';
% folder = 'D:/Dropbox/Dropbox (Vision Lab Cal Tech)/yuli-pharma/';

%% model choices
% method_denoise = 'parametric';
method_denoise = 'svd';
% method_nonlinear = 'relu';
method_nonlinear = 'sigmoid';
[fun_nonlinear, ~] = FitNonlinearity_aet(method_nonlinear);

%% essentials
set(groot,'defaultLineLineWidth',2);

load(['20171207_LNfitparams_' method_denoise '_' method_nonlinear '.mat']);
load(['drugfitparams_' method_denoise '_' method_nonlinear '.mat']);
load('20171207_ethanol_ST.mat');
load('20171207_ethanol_stimulus.mat')
load([folder 'MEA_data/20171207_ethanol/allspk_kilosorted_20171207.mat']);
binsize = 0.1;
fps = 60;
T = 38*29/60*fps;
xq = 1/fps * (1:T);
x = stimulus_barcode;
dx = b(1)*stimulus_barcode+b(2);

%% visualize per cell
for cellnum = selected_cells
    figure('units','normalized','outerposition',[0 0 1 1]); hold on;
    
    switch method_denoise
        case 'parametric'        
            profile_space = paramsspace(:, cellnum);
            profile_time = paramstime(cellnum, :);      
            profile_spacetime = squeeze(paramsl(cellnum, :, :));
            [~, lx] = linear_filter(profile_space, profile_time, x);
            [~, ldx] = linear_filter(profile_space, profile_time, dx);
        case 'svd'
            profile_spacetime = squeeze(paramsl(cellnum, :, :));
            lx = conv_sta(profile_spacetime, x);
            ldx = conv_sta(profile_spacetime, dx);
    end
    pi_nonlinear = paramsnl(cellnum, :);
    output_cntrl = fun_nonlinear(pi_nonlinear, lx);
    output_drug = fun_nonlinear(pi_nonlinear, ldx);

    % control trials
    trials = 5:9;
    [psth_cntrl, edges, rho_cntrl] = PSTH_aet(folder, allspk_bycell, cellnum, trials, binsize, T);    
    psth_cntrl_interp = interpolate_signal(psth_cntrl, edges, xq);

    % drug trials
    trials = 10:22;
    [psth_drug, edges, rho_drug] = PSTH_aet(folder, allspk_bycell, cellnum, trials, binsize, T);    
    psth_drug_interp = interpolate_signal(psth_drug, edges, xq);

    % plot
    scatter(lx, psth_cntrl_interp, 'MarkerFaceColor','b','MarkerEdgeColor','b', 'MarkerFaceAlpha',.2, 'MarkerEdgeAlpha',.2, 'DisplayName', 'control psth');
    scatter(lx, psth_drug_interp, 'MarkerFaceColor','r','MarkerEdgeColor','r', 'MarkerFaceAlpha',.2, 'MarkerEdgeAlpha',.2, 'DisplayName', 'drug psth');
    scatter(ldx, psth_drug_interp, 'MarkerFaceColor','k','MarkerEdgeColor','k', 'MarkerFaceAlpha',.2, 'MarkerEdgeAlpha',.2, 'DisplayName', 'drug model');

    xvals = xlim;
%     thresholds = [pi_nonlinear(2)/pi_nonlinear(1), xvals(1)];
    xdata = xvals(1):0.1:xvals(2);
%     xdata = min(thresholds):0.1:xvals(2);
    output = fun_nonlinear(pi_nonlinear, xdata);
    plot(xdata, output, 'k', 'DisplayName', 'model nonlinearity');
    
    xlabel('model linear filter output');
    ylabel('firing rate (Hz)');
    legend('show');
    legend('Location', 'best');
    
    suptitle(num2str(cellnum));
    
    saveas(gcf, [num2str(cellnum), '_scatter.jpg']);
    close;
end
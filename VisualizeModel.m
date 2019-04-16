folder = '/home/alvita/Dropbox (Vision Lab Cal Tech)/yuli-pharma/';
% folder = 'D:/Dropbox/Dropbox (Vision Lab Cal Tech)/yuli-pharma/';

%% model choices
% method_denoise: parametric, svd
% method_nonlinear: sigmoid, relu
% method_stimfitnl: flicker, barcode
date = '20171207';
method_denoise = 'parametric';
method_nonlinear = 'sigmoid';
method_stimfitnl = '-barcode';
[fun_nonlinear, ~] = nonlinearity_fit(method_nonlinear);

%% load-precomputed
switch date
    case '20171207'
        load([folder 'MEA_data/20171207_ethanol/allspk_kilosorted_20171207.mat']);
    case '20181215'
        load([folder 'MEA_data/20181215_ethanol/20181215_allspkbycell.mat']);
end
folder_results = [folder 'results/' date '/'];
load([folder_results 'ST.mat']);
load([folder_results 'stim.mat']);
load([folder_results 'PSTH.mat']);
load([folder_results 'LN_' method_denoise '-' method_nonlinear method_stimfitnl '.mat']);
load([folder_results 'drug_' method_denoise '-' method_nonlinear method_stimfitnl '.mat']);
folder_figresults = [folder_results 'figs-drug_' method_denoise '-' method_nonlinear method_stimfitnl '/'];
if ~exist(folder_figresults, 'dir')
    mkdir(folder_figresults);
end

fps = 60;
T = 38*29/60;
xq = 1/fps:1/fps:T;
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
    psth_cntrl_interp = interpolate_signal(PSTH.cntrl(cellnum, :), PSTH.edges, xq);
    psth_drug_interp = interpolate_signal(PSTH.drug(cellnum, :), PSTH.edges, xq);
    
    % plot
    scatter(lx, psth_cntrl_interp, 'MarkerFaceColor','b','MarkerEdgeColor','b', 'MarkerFaceAlpha',.2, 'MarkerEdgeAlpha',.2, 'DisplayName', 'control psth');
    scatter(lx, psth_drug_interp, 'MarkerFaceColor','r','MarkerEdgeColor','r', 'MarkerFaceAlpha',.2, 'MarkerEdgeAlpha',.2, 'DisplayName', 'drug psth');
    scatter(ldx, psth_drug_interp, 'MarkerFaceColor','k','MarkerEdgeColor','k', 'MarkerFaceAlpha',.2, 'MarkerEdgeAlpha',.2, 'DisplayName', 'drug model');

    xvals = xlim;
    xdata = xvals(1):0.1:xvals(2);
    output = fun_nonlinear(pi_nonlinear, xdata);
    plot(xdata, output, 'k', 'DisplayName', 'model nonlinearity');
    
    xlabel('model linear filter output');
    ylabel('firing rate (Hz)');
    legend('show');
    legend('Location', 'best');
    
    suptitle(num2str(cellnum));
    
    saveas(gcf, [folder_figresults num2str(cellnum), '.jpg']);
    close;
end
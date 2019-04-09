folder = '/home/alvita/Dropbox (Vision Lab Cal Tech)/yuli-pharma/';
% folder = 'D:/Dropbox/Dropbox (Vision Lab Cal Tech)/yuli-pharma/';

%% model choices
% method_denoise = 'parametric';
method_denoise = 'svd';
% method_nonlinear = 'relu';
method_nonlinear = 'sigmoid';
[fun_nonlinear, ~] = FitNonlinearity_aet(method_nonlinear);

%% essentials
load(['20171207_LNfitparams_' method_denoise '_' method_nonlinear '.mat']);
load('20171207_ethanol_ST.mat');
load('20171207_ethanol_stimulus.mat')
load([folder 'MEA_data/20171207_ethanol/allspk_kilosorted_20171207.mat']);
binsize = 0.1;
fps = 60;
T = 38*29/60*fps;
xq = 1/fps * (1:T);

% % figure settings
set(groot,'defaultLineLineWidth',2)
col = ones(65,3); 
col(1:33,2)=0:1/32:1; 
col(1:33,3)=0:1/32:1;
col(33:65,1) = 1:-1/32:0; 
col(33:65,2) = 1:-1/32:0;
col = flipud(col);

%% select cells to use for drug fitting
THRESH = struct('correlation', 0.5, 'washout', 0.4);
numcells = 153;
all_cells = 1:numcells;
criteria_lnfit = goodnessoffit > THRESH.correlation;
selected_cells = all_cells(criteria_lnfit);
woutcheck = zeros(numcells, 1);
for cellnum = selected_cells
    trials = 5:9;
    [psth_cntrl, ~, rho_cntrl] = PSTH_aet(folder, allspk_bycell, cellnum, trials, binsize, T);    
    trials = 29:38;
    [psth_washout, ~, rho_washout] = PSTH_aet(folder, allspk_bycell, cellnum, trials, binsize, T);    
    woutcheck(cellnum) = mean(psth_washout)/mean(psth_cntrl);
end
criteria_washout = woutcheck > THRESH.washout;
selected_cells = intersect(all_cells(criteria_lnfit), all_cells(criteria_washout));

%% split into cross validation folds
rng(1000);
K = 5;
cvIndices = crossvalind('Kfold', selected_cells, K);

% for k = 1
k = 1;
folderName = ['cvfold' num2str(k) '_' method_denoise '_' method_nonlinear];

if ~exist(folderName, 'dir')
    mkdir(folderName);
end

train_cells = selected_cells(cvIndices ~= k);
test_cells = selected_cells(cvIndices == k);
XY = [];
Z = [];

%% fit single set of transformation parameters for all cells
for cellnum = train_cells
    switch method_denoise
        case 'parametric'        
            profile_space = paramsspace(:, cellnum);
            profile_time = paramstime(cellnum, :);      
            profile_spacetime = squeeze(paramsl(cellnum, :, :));
            [~, lx] = linear_filter(profile_space, profile_time, stimulus_barcode);
        case 'svd'
            profile_spacetime = squeeze(paramsl(cellnum, :, :));
            lx = conv_sta(profile_spacetime, stimulus_barcode);
    end

    % drug trials
    trials = 10:22;
    [psth_drug, edges, rho_drug] = PSTH_aet(folder, allspk_bycell, cellnum, trials, binsize, T);
    psth_drug_interp = interpolate_signal(psth_drug, edges, xq);

    % data for function fitting
    switch method_nonlinear
        case 'relu'
            a1 = paramsnl(cellnum, 1);
            a2 = paramsnl(cellnum, 2);
            m_term = a1 * lx;
            b_term = a1 * sum(profile_spacetime(:)) * ones(size(lx));
            threshold = a2 * ones(size(lx));
            XY = [XY; m_term', b_term', threshold'];
        case 'sigmoid'
            a1 = paramsnl(cellnum, 1);
            a2 = paramsnl(cellnum, 2);
            a3 = paramsnl(cellnum, 3);
            numerator = a1 * ones(size(lx));
            m_term = a2 * lx;
            b_term = a2 * sum(profile_spacetime(:)) * ones(size(lx));
            threshold = a2 * a3 * ones(size(lx));
            XY = [XY; numerator', m_term', b_term', threshold'];
    end
    z = psth_drug_interp;
    Z = [Z; z'];
end

switch method_nonlinear
    case 'relu'
        surfit = @(b, XY) max( b(1)*XY(:, 1) + b(2)*XY(:, 2) - XY(:, 3), 0 );
        b = lsqcurvefit(surfit, [10, 1], XY, Z);
    case 'sigmoid'
        surfit = @(b, XY) XY(1)./( 1+exp( -b(1)*XY(:, 2) - b(2)*XY(:, 3) + XY(:, 4) ) );
        b = lsqcurvefit(surfit, [10, 1], XY, Z);
end
save(['drugfitparams_' method_denoise '_' method_nonlinear '.mat'], 'b', 'selected_cells', 'THRESH');

%% show results of transformation parameters
x = stimulus_barcode;
dx = b(1)*stimulus_barcode+b(2);

for cellnum = selected_cells
    figure('units','normalized','outerposition',[0 0 1 1]);

    subplot(3, 1, 2:3); 
    hold on;

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

    % plot traces
    plot(psth_cntrl_interp, 'b');
    plot(psth_drug_interp, 'r');
    plot(output_cntrl, 'k');
    plot(output_drug, 'k--');
    xlim([0 T]);

    %% cell info
    % STA
    subplot(3, 4, 1);
    sta_vec = ST.average{2}; % kyu version only calc one channel    
    spatial_temporal_map = reshape(sta_vec(:, cellnum)/ max(sta_vec(:, cellnum)), ST.Nx, ST.Nw);
    imagesc(spatial_temporal_map(:, :));
    colormap(col);
    caxis([-1,1]); % set color bar the fixed range of -1 to 1 across trials
    xlabel('Frames')
%     ylabel('Bar position 100um/bar')
    ylabel('Bar Position');
    yticklabels({});

%     if strcmp(method_denoise, 'parametric')
%         subplot(3, 4, 2); hold on;
%         plot(profile_space, 1:numel(profile_space), 'b');
%         plot(profile_space, 1:numel(profile_space), 'k');
%         ylabel('Bar Position');
%         yticklabels({});
%         ylim([1, numel(profile_space)]);
%         set(gca, 'YDir', 'reverse')
% 
%         subplot(3, 4, 3); hold on;
%         plot(0:numel(profile_time)-1, profile_time, 'b');
%         plot(0:numel(profile_time)-1, profile_time, 'k');
%         xlabel('Frames');
%         xlim([0, numel(profile_time)-1]);
%     end

    subplot(3, 4, 4);
    imagesc(profile_spacetime);
    colormap(col);
    caxis([-1, 1]);
    xlabel('Frames')
    %     ylabel('Bar position 100um/bar')
    ylabel('Bar Position');
    yticklabels({});

    suptitle(num2str(cellnum));

    % save figure
    if ismember(cellnum, test_cells) 
        saveas(gcf, [folderName '/test_' num2str(cellnum), '.jpg']);
    else
        saveas(gcf, [folderName '/train_' num2str(cellnum), '.jpg']);
    end
    close;
    disp(cellnum);
end

end
folder = '/home/alvita/Dropbox (Vision Lab Cal Tech)/yuli-pharma/';
% folder = 'D:/Dropbox/Dropbox (Vision Lab Cal Tech)/yuli-pharma/';

%% model choices
% method_denoise: parametric, svd
% method_nonlinear: sigmoid, relu
date = '20171207';
method_denoise = 'parametric';
method_nonlinear = 'sigmoid';
[fun_nonlinear, ~] = nonlinearity_fit(method_nonlinear);

%% essentials
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
load([folder_results 'LN_' method_denoise '-' method_nonlinear '.mat']);

binsize = 0.1;
fps = 60;
T = 38*29/60;
L = T*fps; % 20 seconds at 60 fps
xq = 1/fps * (1:L);

%% figure settings
set(groot,'defaultLineLineWidth',2)
col = ones(65,3); 
col(1:33,2)=0:1/32:1; 
col(1:33,3)=0:1/32:1;
col(33:65,1) = 1:-1/32:0; 
col(33:65,2) = 1:-1/32:0;
col = flipud(col);

%% select cells to use for drug fitting
THRESH = struct('correlation', 0.5, 'washout', 0.4);
numcells = ST.Ncell;
all_cells = 1:numcells;
criteria_lnfit = goodnessoffit > THRESH.correlation;
selected_cells = all_cells(criteria_lnfit);
woutcheck = zeros(numcells, 1);
for cellnum = selected_cells
    woutcheck(cellnum) = mean(PSTH.washout(cellnum, :))/mean(PSTH.cntrl(cellnum, :));
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

% if ~exist(folderName, 'dir')
%     mkdir(folderName);
% end

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

    psth_drug_interp = interpolate_signal(PSTH.drug(cellnum, :), PSTH.edges, xq);    

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
save([folder_results 'drug_' method_denoise '-' method_nonlinear '.mat'], 'b', 'selected_cells', 'THRESH');
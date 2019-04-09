%% get spike data from barcode stimulus
folder = '/home/alvita/Dropbox (Vision Lab Cal Tech)/yuli-pharma/';
% folder = 'D:/Dropbox/Dropbox (Vision Lab Cal Tech)/yuli-pharma/';
date = '20171207';

switch date
    case '20171207'
        load([folder 'MEA_data/20171207_ethanol/allspk_kilosorted_20171207.mat']);
        trial_flckbr = 4;
        trial_barcode = 5;
    case '20181215'
        load([folder 'MEA_data/20181215_ethanol/20181215_allspkbycell.mat']);
end
folder_results = [folder 'results/' date '/'];
load([folder_results 'ST.mat']);
filename = [folder_results 'stim.mat'];

%% reconstruct flicker bar stimulus
stimulus_flckbr = reconstruct_stimulus(trial_flckbr, folder, date);   
  
%% reconstruct barcode stimulus
X = reconstruct_stimulus(trial_barcode, folder, date);
F = size(X, 2);
dx_bc = 4;
 
% match to spatial and temporal resolution of STA
LCR_X = 608;
LCR_Y = 684;
Nx = ST.Nx; % number of X direction checkers
Ny = ST.Ny; % number of Y direction checkers
Nw = ST.Nw; % number of checker map before time zero, determined by the initial w time; 0.4 sec in default 0.4*60hz = 24
fps = 60;
T = 38*29/60*fps; % 20 seconds at 60 fps
stimulus_barcode = zeros(Nx, T);
 
for f = 1:F-1
    % interpolate in space
    stim_bcres = X(:, f);
    stim_bcres = repmat(stim_bcres', dx_bc, 1);
    stim_bcres = stim_bcres(:); % one value per pixel on lcr
 
    increment = LCR_X/Nx;
    stim_stres = zeros(Nx, 1);
    for i = 1:Nx
        left_edge = max(1, round((i-1)*increment));
        right_edge = min(LCR_X, round(i*increment));
        stim_stres(i) = mean(stim_bcres(left_edge:right_edge));
    end
 
    % interpolate in time: 30 frames of each stimulus
    frames = 29;
    left_edge = max(1, (f-1)*frames + 1);
    right_edge = min(T, f*frames);
 
    stimulus_barcode(:, left_edge:right_edge) = repmat(stim_stres, 1, frames);
end
save(filename, 'stimulus_flckbr', 'stimulus_barcode');
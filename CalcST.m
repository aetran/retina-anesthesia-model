folder = '/home/alvita/Dropbox (Vision Lab Cal Tech)/yuli-pharma/';
addpath([folder 'mcd/']);

date = '20171207';
switch date
    case '20171207'
        trial_sta = 4;
        load([folder 'MEA_data/20171207_ethanol/allspk_kilosorted_20171207.mat']); % load spike times from sorted data
        sync_paths = [folder 'MEA_data/20171207_ethanol/20171207_processed/']; % get the sync data of the trial
    case '20181215'
        trial_sta = 3;
        load([folder 'MEA_data/20181215_ethanol/20181215_allspkbycell.mat']);
        sync_paths = [];
end

%% index the corresponding column that is flicker bar trial, check Param to get correct trial
sp_ckbd = allspk_bycell(:, trial_sta);
load([folder 'param_stimuli/' date '/Checkerboard_' num2str(trial_sta) '.mat']); % load stimulus parameters
folder_results = [folder 'results/' date '/'];
if ~exist(folder_results, 'dir')
    mkdir(folder_results);    
end
filename = [folder_results 'ST.mat'];
syncpath1 = [sync_paths 'data' num2str(trial_sta)] ;

sync_time1 = LoadConvertedData(syncpath1);
%debug_synctime  = GetStimTiming(sync_time1);
%debug_synctime  = GetStimTiming_yulimod(sync_time1);
debug_synctime  = GetStimTiming_yulimod_highthres(sync_time1);
figure; (plot(diff(debug_synctime.triggers))); ylim([0 600]); %sanity to check only rep frame is happening no other exception
sync_time1 = GetStimTiming_yulimod(sync_time1); % struct of timing

%% run Hiro's script of STA
stimparam = param; % make sure the param is checkerboard
sync = sync_time1;
disp('done sync time update');

%% run YuLi's STA
clear ST
ST = STanalysisPTBgui_yulimod(sp_ckbd, stimparam, sync);
save(filename, 'ST');

% %% set params
% Nx = ST.Nx; % number of X direction checkers
% Ny = ST.Ny; % number of Y direction checkers
% Nw = ST.Nw; % number of checker map before time zero, determined by the initial w time; 0.4 sec in default 0.4*60hz = 24 frame
% 
% % set colormaps
% col = ones(65,3); 
% col(1:33,2)=0:1/32:1; 
% col(1:33,3)=0:1/32:1;
% col(33:65,1) = 1:-1/32:0; 
% col(33:65,2) = 1:-1/32:0;
% col = flipud(col);
% 
% %% load the STA data, still in vector form
% %RED_STA_VEC = ST.average{1};
% GREEN_STA_VEC = ST.average{2}; % kyu version only calc one channel
% %UV_STA_VEC = ST.average{3};
% disp('done with single channel')
% 
% %% automatically parse every cell and save 
% % make a new dir in save path
% savedir = [folder 'MEA_data/20171207_ethanol/'];
% dir_name = '20171207_STAs';
% if exist(fullfile(savedir, dir_name),'dir')~=7
%     [status,msg] = mkdir(fullfile(savedir, dir_name));
%     assert(status,'Could not create directory ''%s'': %s',fullfile(savedir, dir_name),msg);
% end
% 
% %% reconstrunct the vector STA to XYZ form
% % go thru every trial and save
% 
% for cell_num = 1:size(allspk_bycell,1)
% %for cell_num = 1:5
%     GREEN_STA_MAPs = reshape(GREEN_STA_VEC(:, cell_num)/ max(GREEN_STA_VEC(:, cell_num)), Ny, Nx, Nw); % reshape one cell vector to XY and 3D (69 by 122 by 24)
%     %GREEN_STA_MAPs = reshape(GREEN_STA_VEC(:, cell_num), Ny, Nx, Nw);
%     % also normalized to its max so that ables manual scaling later
% 
%     fig = figure;
%     fig_name = sprintf('20171207 Cell  %d ', cell_num);
% %     suptitle(fig_name);
% 
%     ax_ = subplot(1,1,1);
%     spatial_temporal_map = reshape(GREEN_STA_VEC(:, cell_num)/ max(GREEN_STA_VEC(:, cell_num)), Nx, Nw);
%     imagesc(spatial_temporal_map(:, :)) % chunk out partial
%     colormap(col);
%     caxis([-1,1]); % set color bar the fixed range of -1 to 1 across trials
%     
%     colorbar;
%     xlabel('Frames')
%     ylabel('Bar position 100um/bar')
% 
%     saveas(fig, fullfile(savedir, dir_name, sprintf('%s.png',fig_name)))
%     cla(ax_)
%     close(fig)
% end

function X = reconstruct_stimulus(trial, folder)    

    if nargin == 1
        folder = '/home/alvita/Dropbox (Vision Lab Cal Tech)/yuli-pharma/';
    end

    %% experiment details of all trials
    load([folder 'MEA_data/20171207_ethanol/Param_kilosort_20171207.mat']);
    for i=1:length(Param.DataInfo)
        disp(Param.DataInfo{1, i});
    end

    %% load parameter file for trial
    trial = num2str(trial);

    load([folder 'param_stimuli/20171207/Checkerboard_' trial '.mat'])
    stimparam = param;

    %% get sync data for trial
    sync_paths = [folder 'MEA_data/20171207_ethanol/20171207_processed/'];
    syncpath1 = [sync_paths 'data' trial] ;
    sync_time1 = LoadConvertedData(syncpath1);
    sync_time1 = GetStimTiming_yulimod(sync_time1);
    sync = sync_time1;
    if isfield(sync,'sf')
        sf = sync.sf;
    else
        sf = 10000; % default sampling rate
    end

    %% display parameters
    rng(stimparam.seed, 'twister'); % set the seed
    if ~isfield(stimparam,'rectScreen')
        stimparam.rectScreen = [608 684]; % LCr screen size
    end
    Nx = ceil(stimparam.rectScreen(1)/stimparam.dx);
    Ny = ceil(stimparam.rectScreen(2)/stimparam.dy);
    fps = sf./mean(diff(sync.triggers)); % recorded frame rate in Hz mod in yuli version
    bw = stimparam.bw;
    gauss = stimparam.gauss;

    t = sync.triggers/sf ; % trigger timings in sec; in Yuli's version, trigger is precalculated in getstimTiming
    N = length(t)-1; % total number of frames

    %% stimulus reconstruction
    j = round(linspace(0,N,11));
    j = j(2:end-1);
    tic;

    X = zeros(Nx, N);

    for i=1:N
    % for i = 1
        %% create noise image frame: see also Checkerboard.m for PTBgui
        noiseimg = randn(Nx*Ny,3);
        if gauss > 0  % gaussian
            noiseimg = gauss*noiseimg; % + gray
        else % binary
            for k=1:3 
                noiseimg(noiseimg(:,k)<=0,k)=-1; % black
                noiseimg(noiseimg(:,k)~=-1,k)=1; % white
            end
        end
        noiseimg = reshape(noiseimg,Ny,Nx,3);
        if bw==1
            for k=2:3
                noiseimg(:,:,k) = noiseimg(:,:,1);
            end % black-and-white
        elseif any(bw==-3:-1)
            for k=setxor(abs(bw),1:3)
                noiseimg(:,:,k)=0;
            end % single-color
        end 

        X(:, i) = noiseimg(:, :, 1);

        %% progress indicator
        if any(i==j)
            fprintf('Computing %d out of %d frames. ', i, N);
            toc;
        end
    end
end
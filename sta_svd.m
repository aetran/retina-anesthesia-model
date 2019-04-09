function [profile_spacetime] = sta_svd(cellnum, ST, ntimepoints_to_keep, npixels_to_keep, components_to_use, ncomponents_to_disp, opt_plot)
if nargin == 5
    opt_plot = 0;
end

ntimepoints_to_keep = min(ntimepoints_to_keep, ST.Nw);

% red blue colormap
col = ones(65,3); 
col(1:33,2)=0:1/32:1; 
col(1:33,3)=0:1/32:1;
col(33:65,1) = 1:-1/32:0; 
col(33:65,2) = 1:-1/32:0;
col = flipud(col);

% cell's STA
sta_full = ST.average{2}(:, cellnum);
sta_full = reshape(sta_full, ST.Ny, ST.Nx, []);

% crop profile in time
sta_crop = sta_full(:, :, end-ntimepoints_to_keep+1:end);

% max pixel
sta_crop = reshape(sta_crop, ST.Ny*ST.Nx, []) / max(abs(sta_crop(:)));
[max_pos,max_time]=find(abs(sta_crop)==1); % max_time gives the max of temporal profile
if length(max_time)>1
    max_time=max_time(1);
end
a=mod(max_pos, ST.Ny);
if a==0
    a=ST.Ny;
end
b=(max_pos-a)/ST.Ny+1; 
a=a(1);
b=b(1);
sta_crop = reshape(sta_crop,ST.Ny,ST.Nx,[]);

% crop profile in space
y_last=min(ST.Ny,a+npixels_to_keep);
y_first=max(1,a-npixels_to_keep);
x_last=min(ST.Nx,b+npixels_to_keep);
x_first=max(1,b-npixels_to_keep);
new_Ny=y_last-y_first+1;
new_Nx=x_last-x_first+1;
sta_crop=sta_crop(y_first:y_last, x_first:x_last, :);

% figure
if opt_plot 
    figure; 
    subplot(2, ncomponents_to_disp+1, 1); % raw: spatial profile at max pixel
    imagesc(squeeze(sta_crop(:,:,max_time)), [-1 1]);
    colormap(col);
    axis equal tight off
    title('Raw');

    subplot(2, ncomponents_to_disp+1, ncomponents_to_disp+2); % raw: time course of max pixel
    tv=linspace(0-(ntimepoints_to_keep-1)*1/ST.fps,0,ntimepoints_to_keep);
    plot(tv, squeeze(sta_full(a, b, end-ntimepoints_to_keep+1:end)));
    xlim([tv(1),tv(end)]);
    xlabel('Time to spike (s)');
    box off
end

%% compute SVD
sta_crop = reshape(sta_crop,new_Ny*new_Nx,[]);
[u,s,v]=svd(sta_crop,'econ');

if opt_plot % coded for 1D stimulus
    for i = 1:ncomponents_to_disp
        z=reshape(u(:,i), new_Ny, new_Nx, []);
    %     zmax = max(abs(u(:,i)));
        zmax = max(z(:));

        subplot(2, ncomponents_to_disp+1, i+1);
        imagesc(z, zmax*[-1 1]);
        colormap(col);
        axis equal tight off
        title(['SVD component ' num2str(i)]);

        subplot(2, ncomponents_to_disp+1, ncomponents_to_disp+i+2);
        plot(tv, v(:,i));
        xlabel('Time to spike (s)');
        xlim([tv(1), tv(end)]);
        box off
    end
    suptitle(num2str(cellnum));
    set(gcf, 'Position', get(0, 'Screensize'));
    filename = [num2str(cellnum) '_svd_details'];
    print(gcf, filename, '-dpng');
end

%% smoothed spatiotemporal receptive field
profile_crop = u(:, components_to_use) * s(components_to_use, components_to_use) * v(:, components_to_use)';
profile_crop = reshape(profile_crop, new_Ny, new_Nx, ntimepoints_to_keep);
profile_full = zeros(ST.Ny, ST.Nx, ST.Nw);
profile_full(y_first:y_last, x_first:x_last, end-ntimepoints_to_keep+1:end) = profile_crop;
profile_spacetime = squeeze(profile_full);

% % match to original magnitude
% scale = 1/max(abs(profile_full(:)));
% profile_spacetime = profile_spacetime*scale;

% figure, coded for 1D
if opt_plot
    figure;
    
    % each component
    for j = 1:ncomponents_to_disp
        z = u(:, j) * s(j, j) * v(:, j)';
        z = reshape(z, new_Ny, new_Nx, ntimepoints_to_keep);
        z = squeeze(z);
        
        subplot(2, ncomponents_to_disp, j);
        imagesc(z, [-1 1] );
        axis equal
        box on
        set(gca, 'TickLength',[0 0])
        set(gca, 'XTickLabel',[])
        set(gca, 'YTickLabel',[])
        colormap(col);
        title(num2str(j));
    end
    
    % raw STA
    subplot(2, ncomponents_to_disp, j+1); 
    sta_crop = reshape(sta_crop, new_Ny, new_Nx, ntimepoints_to_keep);
    sta_crop = squeeze(sta_crop);
    imagesc(sta_crop, [-1 1]);
    axis equal
    box on
    set(gca, 'TickLength',[0 0])
    set(gca, 'XTickLabel',[])
    set(gca, 'YTickLabel',[])
    colormap(col);
    title('Raw');
    
    subplot(2, ncomponents_to_disp, j+2);
    sta_full = reshape(sta_full, ST.Ny, ST.Nx, ST.Nw);
    sta_full = squeeze(sta_full);
    imagesc(sta_full, [-1 1] );
    axis equal
    box on
    set(gca, 'TickLength',[0 0])
    set(gca, 'XTickLabel',[])
    set(gca, 'YTickLabel',[])
    colormap(col);
    title('Raw (Full)');
   
    subplot(2, ncomponents_to_disp, j+3);
    profile_full = reshape(profile_full, ST.Ny, ST.Nx, ST.Nw);
    profile_full = squeeze(profile_full);
    imagesc(profile_full, [-1 1]);
    axis equal
    box on
    set(gca, 'TickLength',[0 0])
    set(gca, 'XTickLabel',[])
    set(gca, 'YTickLabel',[])
    colormap(col);
    title('Denoised (Full)');
    
    suptitle(num2str(cellnum));
    set(gcf, 'Position', get(0, 'Screensize'));
    filename = [num2str(cellnum) '_svd_frames'];
    print(gcf, filename, '-dpng');
end

end
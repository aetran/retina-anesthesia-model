function plot_raster(spikes, row, color)
    if nargin == 1
        row = 1;
    end
    if nargin <= 2
        color = 'b';
    end
    hold on;
    
    if ~isempty(spikes)
        spikes = [spikes(1); spikes];
        y = row*ones(size(spikes));
        y(1:2:end) = row-1;
        stairs(spikes, y, color); % spikes
        plot([min(spikes) max(spikes)], [row-1 row-1], 'w', 'LineWidth', 2); % erase lines between rows
        plot([min(spikes) max(spikes)], [row row], 'w', 'LineWidth', 2);
    end
end
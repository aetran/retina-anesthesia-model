function [raw_space, raw_time, profile_space, profile_time, profile_spacetime] = sta_parametric(cellnum, ST)

    %% get STA of specific cell
    sta_vec = ST.average{2}(:, cellnum);
    sta = reshape(sta_vec / max(abs(sta_vec)), ST.Nx, ST.Nw); % Nx x Nw, omit Ny because one dimensional

    %% find spatial and temporal profiles at max pixel
    [~, ind] = max(abs(sta(:)));
    [i, j] = ind2sub(size(sta), ind);
    raw_space = sta(:, j);
    if sta(ind) < 0 
        raw_space = -raw_space;
    end
    raw_time = fliplr(sta(i, :));
    
    %% smooth
    xdata = 1:numel(raw_space);
    ydata = raw_space';
    fun_gaussian = @(x, xdata) (x(1) * exp(-((xdata - x(2)).^2 ) / 2*(x(3))^2 ) );
    [~, peak] = max(raw_space);
    x0 = [1, peak, 1];
    pi_gaussian = lsqcurvefit(fun_gaussian, x0, xdata, ydata);
    profile_space = fun_gaussian(pi_gaussian, xdata);
    profile_space = profile_space';
    
    xdata = (0:23)*1000/60;
    ydata = raw_time;
    fun_temporal = @(x, xdata) ((xdata/x(1)).^x(2)).*exp(-x(2)*(xdata/x(1) -1))...
        -x(3)*((xdata/x(4)).^x(5)).*exp(-x(5)*(xdata/x(4) -1));
    if max(raw_time) > abs(min(raw_time))
        x0 = [30, 5, 1, 100, 10];
    else
        x0 = [60, 20, 10, 50, 30];
    end
    pi_temporal = lsqcurvefit(fun_temporal, x0, xdata, ydata);
    profile_time = fun_temporal(pi_temporal, xdata);
    
    profile_spacetime = fliplr(profile_time) .* profile_space;
    
end
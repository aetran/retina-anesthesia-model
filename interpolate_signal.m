function psth_new = interpolate_signal(psth, edges, edges_new)

x = edges(1:end-1);
xq = edges_new;
vq = interp1(x, psth, xq);
psth_new = vq;
psth_new(isnan(psth_new)) = 0;

end
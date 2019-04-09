function output_time = conv_sta(sta, stimulus)
    T = size(stimulus, 2);
    f = size(sta, 2);    
    output_time = zeros(1, T);
    for t = f:T
        output_time(t) = sum(sum(sta.*stimulus(:, t-f+1:t)));
    end
end
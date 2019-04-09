function [output_space, output_time] = linear_filter(profile_space, profile_time, stimulus)

output_space = profile_space' * stimulus;
output_time = conv(output_space, profile_time);
output_time = output_time(1:size(stimulus, 2));
% output_time = output_time(size(profile_time, 2):end);

end
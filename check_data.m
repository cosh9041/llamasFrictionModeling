function [faults, indicies] = check_data( data_x, data_y, step, N_ma, threshold )
% Checks data for faults, checking every 'step' steps, and calculating
% fault based on last 'N_ma' number of moving average points. Triggers
% fault if above a certain threshold. 
% Assumes data_x and data_y have same length. 

% initialize vars
faults = zeros(floor(length(data_x)/step) - 1, 1);
indicies = zeros(floor(length(data_x)/step) - 1, 1);

j = 1;
for i = 1:step:(length(data_x) - step)
    faults(j) = check_threshold_residual(data_x(1:i), ...
        data_y(1:i), N_ma, threshold);
    indicies(j) = i;
    j = j+1;
end
end


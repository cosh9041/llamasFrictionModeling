function is_fault = check_threshold_residual( data_x, data_y, sz_ma, threshold )
is_fault = false;
start_x = length(data_x) - sz_ma;
start_y = length(data_y) - sz_ma;
if start_x < 1
    start_x = 1;
end
if start_y < 1
    start_y = 1;
end
delta = (data_x(start_x:end) - data_y(start_y:end));
mean_delta = abs(mean(delta));
std_delta = std(delta);

if mean_delta > threshold
    if (mean_delta - threshold) > std_delta
        is_fault = true;
        return;
    end
end
end


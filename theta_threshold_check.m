function theta_threshold_check( )
%Creates mock faulty sensor data and runs it through threshold fault
%checking algorithm

%% theta shiz
N = 1000;
d_theta_coarse = (rand(N, 1) - 0.5)*2 + 0.1*randn(N,1) + ...
    [0*ones(N/2, 1); 5*ones(N/2, 1)]; % 5 degree static offset, manifested in 
d_theta_fine = (rand(N, 1) - 0.5)*1 + 0.1*randn(N,1); % d_theta reading
t = linspace(0, N/50, N);
fig = figure
hold on
step = 10;
[fault_status, indicies] = check_data(d_theta_coarse, d_theta_fine, step, ...
    50, (std(d_theta_fine) + std(d_theta_coarse)));

fault_line = find(fault_status == 1);
fault_start = indicies(fault_line(1));
fault_end = indicies(end);


dim_bad = [t(fault_start)/t(end), .11, ...
    (t(fault_end) - t(fault_start))/t(end)-.075,1-.180];
dim_good = [0.13, 0.11, t(fault_start)/t(end) - .13,1-.180];
annotation('rectangle',dim_bad,'FaceColor','r','FaceAlpha',.2, 'Units', 'pixels')
annotation('rectangle',dim_good,'FaceColor','g','FaceAlpha',.2)
% annotation('textbox',[0, ,'String',str,'FitBoxToText','on');

plot(t, d_theta_coarse)
plot(t, d_theta_fine)
plot(t(indicies), fault_status)
end

% This commented out code will plot the locations of moving average for
% example purposes
% j = 1;
% for i = 1:step:(length(d_theta_coarse) - step)
% %     idx_end = (i-1) + step;
% %     ma_coarse(j) = mean(d_theta_coarse(i:idx_end));
% %     std_coarse(j) = std(d_theta_coarse(i:idx_end));
% %     t_ma(j) = mean(t(i:idx_end));
% %     ma_fine(j) = mean(d_theta_fine(i:idx_end));
% %     std_fine(j) = std(d_theta_fine(i:idx_end));
%     faults(j,1) = check_threshold_residual(d_theta_coarse(1:i), ...
%         d_theta_fine(1:i ), 50, (std(d_theta_fine) + std(d_theta_coarse)));
%     locations(j) = t(i);
%     j = j+1;
% end
% errorbar(t_ma, ma_coarse, std_coarse, 'x--')
% errorbar(t_ma, ma_fine, std_fine, 'x--')
% plot(t(indicies), fault_status*1)
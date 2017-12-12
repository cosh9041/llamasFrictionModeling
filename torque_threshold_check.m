function [fault_status, indicies] = ...
    torque_threshold_check( f, f_hat, threshold, res_std, omegas )
% Creates mock friction torque data and runs it through fault check
% threshold alrogithm. Mock data friction function is based off of ASEN
% 3200 spin module data. 
N = 1000;
if nargin < 5
    omegas = 10*rand(N,1);
end
t =  linspace(0, 10, N);
N = length(omegas);
nominal = f(1, omegas) + randn(length(omegas), 1)*res_std;
induced = [f(1, omegas(1:N/2)) + randn(N/2, 1)*res_std; ...
    f_hat(1, omegas(N/2+1:end)) + randn(N/2, 1)*res_std];

step = 10;
[fault_status, indicies] = check_data(nominal, induced, step, 50, threshold);

figure
plot(t, nominal, '*')
hold on
plot(t, induced, 'x')

fault_line = indicies(find(fault_status == 1));
if isempty(fault_line)
    fault_start = N;
    dim_bad = [0 0 0 0];
    dim_good = [0.13, 0.11, t(fault_start)/t(end) - .13 - .095,1-.185];
else
    fault_start = fault_line(1);
    dim_bad = [t(fault_start)/t(end), .11, ...
        (t(end) - t(fault_start))/t(end)-.095,1-.185];
    dim_good = [0.13, 0.11, t(fault_start)/t(end) - .13,1-.185];
end

annotation('rectangle',dim_bad,'FaceColor','r','FaceAlpha',.2)
annotation('rectangle',dim_good,'FaceColor','g','FaceAlpha',.2)
plot(t(indicies), fault_status)
end


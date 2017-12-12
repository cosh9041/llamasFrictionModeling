function [fault_status, indicies] = ...
    torque_threshold_check( f, f_hat, threshold, res_std, omegas )
% Creates mock friction torque data and runs it through fault check
% threshold alrogithm. Mock data friction function is based off of ASEN
% 3200 spin module data. 
if nargin < 5
    omegas = 10*rand(1000,1);
end
nominal = f(1, omegas) + randn(length(omegas), 1)*res_std;
induced = f_hat(1, omegas) + randn(length(omegas), 1)*res_std;

step = 10;
[fault_status, indicies] = check_data(nominal, induced, step, 50, threshold);
end


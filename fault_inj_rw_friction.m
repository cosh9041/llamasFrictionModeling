%% Modeling for FI PDR in reaction wheels
clear all; close all; clc;

tspan = [0, 20];
I = 0.01; 


[f, p, omega_nom, t_f_nom] = friction_func(I);

f_hat = @(t, omega) f(t, omega) + 5.9*p(2)*sign(omega); % INcreasing friction due to bearing failure by 50%
f_threshold = @(t, omega) f(t, omega) + 4*p(2)*sign(omega);

omega0 = [0, 1, 10:10:160];
alpha_diff_omega = zeros(1, length(omega0));
t_adtl_max = zeros(1, length(omega0));
t_c = 2;
t_induced = zeros(1, length(omega0));
for i = 1:length(omega0) 
    [alpha, omega, t_induced] = friction_alpha(f, f_hat, I, t_c, omega0(i), tspan);
    t_induced0(i) = t_induced(1);
end

%% calculate normal data deviation from nominal frition function
residual = f(1, omega_nom) - t_f_nom;
res_std = std(residual);
E = ones(length(residual), 1).*res_std;
figure
hold on
t_induced = interp1(omega0, t_induced0, omega_nom) + randn([length(omega_nom), 1])*res_std;
errorbar(omega_nom, t_induced, E, '*--', 'LineWidth', 1.5)
plot(omega0, f_threshold(1, omega0), '--', 'LineWidth', 1.5);
% plot(omega0, f(1, omega0), 'LineWidth', 1.5)
errorbar(omega_nom, t_f_nom, E, 'p--');
legend({'\tau_{f,measured}', '\tau_{f,threshold}', '\tau_{f,nominal}', 'tau_{f,act}'},'Location', 'east', 'FontSize', 16)
xlabel('\omega (rad/s)', 'FontSize', 16)
ylabel('\tau_f (N-m)', 'FontSize', 16)
title('Comparing sensed friction torque vs nominal and threshold', 'FontSize', 16)


%% theta shiz
N = 1000;
d_theta_coarse = (rand(N, 1) - 0.5)*5 + 0.4*randn(N,1) + 5; % 5 degree static offset, manifested in 
d_theta_fine = (rand(N, 1) - 0.5)*2 + 0.2*randn(N,1); % d_theta reading
t = linspace(1, N/50, N);
figure
hold on
plot(t, d_theta_coarse)
plot(t, d_theta_fine)
step = 50;
j = 1;
for i = 1:step:(length(d_theta_coarse) - step)
    idx_end = (i-1) + step;
    ma_coarse(j) = mean(d_theta_coarse(i:idx_end));
    std_coarse(j) = std(d_theta_coarse(i:idx_end));
    t_ma(j) = mean(t(i:idx_end));
    ma_fine(j) = mean(d_theta_fine(i:idx_end));
    std_fine(j) = std(d_theta_fine(i:idx_end));
    j = j+1;
end
errorbar(t_ma, ma_coarse, 2*std_coarse, 'x--')
errorbar(t_ma, ma_fine, 2*std_fine, 'x--')

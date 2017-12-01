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
%     figure
%     hold on
%     plot(alpha.time_induced, t_induced, '*')
%     plot(omega.t_induced, f_threshold(1, omega.induced));
% %     tStr = ['\Delta\alpha vs time at \omega = ', sprintf('%.2f rad/s', omega0(i))];
% %     title(tStr)
%     xlabel('Time (s)')
%     ylabel('\tau_f (N-m))')
%     legend('\hat{tau}_f', '\tau_{threshold}')
% %     t_adtl_max(i) = max(t_adtl);
%     t_max_loc(i) = find(max(t_adtl));
    t_induced0(i) = t_induced(1);
end

%% calculate normal data deviation from nominal frition function
residual = f(1, omega_nom) - t_f_nom;
res_std = std(residual);
figure
hold on
t_induced = interp1(omega0, t_induced0, omega_nom) + randn([length(omega_nom), 1])*res_std;
plot(omega_nom, t_induced, '*--', 'LineWidth', 1.5)
plot(omega0, f_threshold(1, omega0), '--', 'LineWidth', 1.5);
plot(omega0, f(1, omega0), 'LineWidth', 1.5)
plot(omega_nom, t_f_nom, 'p--');
% plot(omega_nom, t_induced, '*--');
% plot(omega_nom, t_f_nom)
% axis([-omega, omega0(end), 0, max(t_induced) + 0.2*max(t_induced)]);
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
% plot(omega0, alpha_diff_omega, '*')
% 
% t_c = 0:0.25:5;
% omega0 = 0;
% alpha_diff_tc = zeros(1, length(t_c));
% 
% for j = 1:length(t_c) 
%     [alpha, omega, t_adtl] = friction_alpha(f, f_hat, I, t_c(j), omega0, tspan);
%     figure
%     hold on
%     plot(alpha.t_act, t_adtl, '*')
%     plot([0, alpha.t_act(end)], [4*p(2), 4*p(2)]);
%     tStr = ['\tau_{adtl} vs time at \tau_c = ', sprintf('%.2f', t_c(j))];
%     title(tStr)
%     xlabel('Time (s)')
%     ylabel('\tau_{adtl}')
%     legend('\tau_{adtl}', '\tau_{threshold}')
%     t_adtl_max(j) = max(t_adtl);
% end
% figure
% hold on
% plot(t_c, t_adtl_max, '*')
% plot([0, t_c(end)], [4*p(2), 4*p(2)], '--', 'LineWidth', 1.5);
% axis([0, t_c(end), 0, max(t_adtl_max) + 0.2*max(t_adtl_max)]);
% legend('\tau_{adtl}', '\tau_{threshold}')
% xlabel('\tau_c [N-m]')
% ylabel('\tau_{adtl}')
% title('Added friction torque while holding \omega_0 = 0')

% figure
% hold on
% plot(t_c, alpha_diff_tc, '*')
% title('\Delta \alpha vs \tau_c')
% xlabel('\tau_c (N-m)')
% ylabel('\Delta \alpha (rad/s^2)')

% 
% figure
% hold on
% plot(omega, f(t, omega));
% plot(omega, f_hat(t, omega));
% legend('\tau_f nominal', '\tau_f induced')
% title('Induced vs nominal friction')
% 
% figure
% hold on
% plot(t, omega)
% plot(t_hat, omega_hat)
% legend('Omega w/ nominal friction', 'Omega w/ induced friction')
% xlabel('Time (s)')
% ylabel('Reaction wheel speed (rad/s)')
% tStr = ['Reaction wheel speed vs time for \tau_c = ', sprintf('%.2f N-m', t_c)];
% title(tStr)
% 
% figure
% hold on
% plot(alpha.t_act, alpha.act);
% plot(alpha.t_induced, alpha.induced)
% legend('\alpha w/ nominal friction', '\alpha w/ induced friction', ...
%        '\alpha w/ nominal friction, eval with data', '\alpha w/ induced friction, eval w/ data')
% xlabel('Time (s)')
% ylabel('\alpha (rad/s^2)')
% tStr = ['\alpha vs time for \tau_c = ', sprintf('%.2f N-m', t_c)];
% title(tStr)
function print_torque_plot( f, f_hat, I, tspan, omega_nom, t_f_nom )

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
% plot(omega0, f_threshold(1, omega0), '--', 'LineWidth', 1.5);
plot(omega0, f(1, omega0), 'LineWidth', 1.5)
errorbar(omega_nom, t_f_nom, E, 'p--');
legend({'\tau_{f,induced}', ...%'\tau_{f,threshold}', 
    'f(\omega)', '\tau_{f,actual}'},'Location', 'east', 'FontSize', 16)
xlabel('\omega (rad/s)', 'FontSize', 16)
ylabel('\tau_f (N-m)', 'FontSize', 16)
title('Comparing sensed friction torque vs nominal and threshold', 'FontSize', 16)
end


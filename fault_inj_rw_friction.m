%% Modeling for FI PDR in reaction wheels
clear all; close all; clc;

tspan = [0, 20];
I = 0.01; %Estimated MOI of ASEN 3200 reaction wheels, kg-m^2

[f, p, omega_nom, t_f_nom] = friction_func(I);

f_hat = @(t, omega) f(t, omega) + 6*p(2)*sign(omega); % INcreasing friction due to bearing failure by 50%
f_threshold = @(t, omega) f(t, omega) + 4*p(2)*sign(omega);
residual_threshold = 4*p(2);
print_torque_plot(f, f_hat, I, tspan, omega_nom, t_f_nom);

residual = f(1, omega_nom) - t_f_nom;
res_std = std(residual);
E = ones(length(residual), 1).*res_std;

%% Run threshld check on mock theta data, given 5 degree offset
theta_threshold_check();

%% Run threshold check on induced friction vs nominal
torque_threshold_check(f, f_hat, residual_threshold, res_std);
%% Run threshold check on nominal vs nominal friction (Should not fault)
torque_threshold_check(f, f, residual_threshold, res_std);

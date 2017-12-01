function [ alpha, omega, t_induced ] = friction_alpha( nat_f, induced_f, I, t_c, omega0, tspan)
%
w_dot = @(t, omega) (t_c - nat_f(t, omega))/I;
w_dot_hat = @(t, omega) (t_c >= induced_f(t, omega))*(t_c - induced_f(t, omega))/I;

[omega.t_act, omega.act] = ode45(w_dot, tspan, omega0);
[omega.t_induced, omega.induced] = ode45(w_dot_hat, tspan, omega0);

alpha.act = diff(omega.act)./diff(omega.t_act);
alpha.induced = diff(omega.induced)./diff(omega.t_induced);
alpha.t_act = omega.t_act(1:end-1) + .5 * diff(omega.t_act);
alpha.time_induced = omega.t_induced(1:end-1) + .5 * diff(omega.t_induced);

% t_cdf = I*(alpha.act - interp1(alpha.t_induced, alpha.induced, alpha.t_act));
t_induced = t_c - I*(alpha.induced);
end


function [nat_f, p, omega, t_f] =  friction_func( I )
%Calculates the natural friction in reaction wheel data.
for i = 5
    filestr = sprintf('AM_020317_18A_RWHEEL_COMP_0%.0f', i);
    rwheel(:,:,i) = load(filestr);
    rwheel_out(:,:,i) = rwheel(rwheel(:,2,i) ~= 0, :, i);
    alpha = diff(rwheel_out(:,3,i))./diff(rwheel_out(:,1,i));
    t_meas = I*alpha;
    t_meas_time = rwheel_out(1:end-1,1,i) + 0.5*diff(rwheel_out(:,1,i));
    t_diff(:,i) = interp1(rwheel_out(:,1,i), rwheel_out(:,2,i), t_meas_time) - t_meas;
    rwheel_diff(:,i) = interp1(rwheel_out(:,1,i), rwheel_out(:,3,i), t_meas_time);
    omega = rwheel_diff(:,i);
    t_f = t_meas;
    
    figure
    hold on
    plot(rwheel_out(:,3,i), rwheel_out(:,2,i));
    plot(omega, t_meas, '--*');
    tStr = ['Torque profile of ASEN 3200 Spin Modules'];
    title(tStr, 'FontSize', 16)
    ylabel('\tau (N-m)', 'FontSize', 16)
    xlabel('\omega (rad/s)', 'FontSize', 16)

    t_diff(t_diff(:,i) < 0, i) = 0;
    plot(omega, t_diff(:,i), 'k--+')
    [p] = polyfit(rwheel_diff(:,i), t_diff(:,i), 1);
    nat_f = @(t, omega) p(1)*omega + p(2)*sign(omega);
    plot(omega, nat_f(t_meas_time, rwheel_diff(:,i)), 'k', 'LineWidth', 1.5)
    legend({'\tau_c', '\tau_{plant}', '\tau_f', 'Fitted \tau_f line'})
%     figure
%     hold on
%     plot(omega, t_diff(:,i), '*')
%     plot(omega, nat_f(t_meas_time, rwheel_diff(:,i)))
%     f_hat = @(t, omega) nat_f(t, omega) + 5.9*p(2);
%     plot(rwheel_diff(:,i), f_hat(t_meas_time, rwheel_diff(:,i)))
%     title('Nominal friction torque in ASEN 3200 spin module reaction wheels')
%     legend('Raw Data', 'Curve fit', 'Failure Friction')
%     xlabel('\omega (rad/s)')
%     ylabel('\tau_f (N-m)')
    

end

end


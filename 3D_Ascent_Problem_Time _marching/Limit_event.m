function [value, isterminal, direction] = Limit_event(t, x)
    % Define event function to stop integration at 400 km altitude
    Re = 6371000;                   % Radius of earth in meters
    value = x(1) - (Re + 400000);    % 400 km altitude
    isterminal = 1;  % Stop the integration
    direction = 0;  % Any direction
end
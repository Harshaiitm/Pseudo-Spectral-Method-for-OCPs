function velocity_guess = rocket_trajectory_guess(lat, lon, alt_i, v_i, alt_f, inc_f, n_points)
    % Constants
    a = 6378137.0; % WGS-84 semi-major axis in meters
    e = 0; % WGS-84 first eccentricity
    omega = 7.2921150e-5; % Earth's angular velocity in rad/s
    R_earth = 6371000; % Earth's mean radius in meters

    % Convert lat, lon to radians
    lat = deg2rad(lat);
    lon = deg2rad(lon);
    inc_f = deg2rad(inc_f);

    % Initial position in ECEF
    N = a / sqrt(1 - e^2 * sin(lat)^2);
    X_i = (N + alt_i) * cos(lat) * cos(lon);
    Y_i = (N + alt_i) * cos(lat) * sin(lon);
    Z_i = ((1 - e^2) * N + alt_i) * sin(lat);

    % Final altitude in meters
    alt_f = alt_f * 1000;

    % Final position in ECEF (assuming same lat, lon for simplicity)
    N_f = a / sqrt(1 - e^2 * sin(lat)^2);
    X_f = (N_f + alt_f) * cos(lat) * cos(lon);
    Y_f = (N_f + alt_f) * cos(lat) * sin(lon);
    Z_f = ((1 - e^2) * N_f + alt_f) * sin(lat);

    % Linspace to generate initial guesses
    v_x_guess = linspace(0, 10 * v_i, n_points);
    v_y_guess = linspace(0, 10 * v_i, n_points);
    v_z_guess = linspace(0, 10 * v_i, n_points);

    % Preallocate velocity guess array
    velocity_guess = zeros(n_points, 3);

    % Loop through and fill velocity guesses
    for i = 1:n_points
        % Initial ENU velocities (assuming initial velocity is purely vertical for simplicity)
        v_E = v_x_guess(i);
        v_N = v_y_guess(i);
        v_U = v_z_guess(i);

        % Transformation matrix from ENU to ECEF
        T = [-sin(lon), -sin(lat)*cos(lon), cos(lat)*cos(lon);
              cos(lon), -sin(lat)*sin(lon), cos(lat)*sin(lon);
              0,         cos(lat),          sin(lat)];

        % Transform ENU velocity to ECEF velocity
        v_ecef_enu = [v_E; v_N; v_U];
        v_ecef_transformed = T * v_ecef_enu;

        % Velocity due to Earth's rotation
        v_rotation = omega * [-Y_i; X_i; 0];

        % Final ECEF velocity
        v_ecef = v_ecef_transformed + v_rotation;

        % Store in array
        velocity_guess(i, :) = v_ecef';
    end
end



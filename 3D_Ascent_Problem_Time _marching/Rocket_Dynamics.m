function dxdt = Rocket_Dynamics(t, x)
    r = x(1);    % Radial position   
    theta = x(2);   % Longitude
    phi = x(3);   % Latitude     
    V = x(4);     % Velocity
    gamma = x(5);   % Flight path angle
    psi = x(6);     % Azimuth
     m = x(7);

    % Density Model Call
    Re = 6371000;                   % Radius of earth in meters
    h_scale = 8500;             
    mu = 3.986012e14;               % Gravitational parameter "GM" in m^3/s^2
    Omega = 2*pi/(24*60*60);      % Sideral Rotation Rate (rad/s)
    rho0 = 1.225;                   % air density at Sea level 
    g0 = 9.80665;                   % acceleration due to gravity at sea level
    m0_1= 248950;                   % 1st stage total mass
    m0_2= 134040;                   % 2nd stage total mass
    m0 = m0_1 + m0_2;
    %m = m0;

    g = mu / r^2;
    rho = rho0 * exp(-(r - Re) / h_scale);
    Cd = 0.3;
    Cl = 0.5;
    A_ref = 61;
    Thrust = m * g * 1.2;
    Isp = 300;

    D = 0.5 * rho * V^2 * A_ref * Cd;
   
    L = 1/2 * rho * V^2 * A_ref * Cl;

    sigma = 0;

    Ft = Thrust - D;
    Fn = L;

    dxdt = [
         V * sin(gamma); 
         V * cos(gamma) * cos(psi) / (r * cos(phi));
         V * cos(gamma) * sin(psi) / r; 
           Ft / m - g * sin(gamma) + Omega^2 * r * cos(phi) * (sin(gamma) * cos(phi) - cos(gamma) * sin(psi) * sin(phi));
        ((Fn / m) * cos(sigma) - g * cos(gamma) + (V^2 * cos(gamma)) / r + 2 * Omega * V * cos(psi) * cos(phi) + Omega^2 * r * cos(phi) * (cos(gamma) * cos(phi) + sin(gamma) * sin(psi) * sin(phi))) / V;
        (((Fn / m) * (sin(sigma) / cos(gamma))) - (V^2 / r) * cos(gamma) * cos(psi) * tan(phi) + 2 * Omega * V * (tan(gamma) * sin(psi) * cos(phi) - sin(phi)) - (Omega^2 * r * cos(psi) * sin(phi) * cos(phi)) / cos(gamma)) / V;
        -Thrust/(g0*Isp)
    ];
end
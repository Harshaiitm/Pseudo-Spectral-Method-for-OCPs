function [R1_I,V1_I,qn1] = lat_long_elev_azi_vec1(M,problem) 

Elev_i = problem.Elev_i;
Azim_i = problem.Azim_i;
hi = problem.hi;
Vi = problem.Vi;
Elev_s = problem.Elev_s;
Azim_s = problem.Azim_s;
hf_s = problem.hf_s;
Vf_s = problem.Vf_s;
Re = problem.Re;
nodes = problem.nodes;
D = problem.D;

t0 = 0;
stage_time = 160;

load time_alt_lat_long.csv
time = time_alt_lat_long(1:9,1);
alt = time_alt_lat_long(1:9,2)*1000;
lat = time_alt_lat_long(1:9,3);
long = time_alt_lat_long(1:9,4);

t_1 = ((stage_time-t0)/2).*nodes+(stage_time+t0)/2;

Altitude0_1 = interp1(time, alt, t_1, 'pchip');
latitude_1 = deg2rad(interp1(alt, lat, Altitude0_1, 'pchip'));
longitude_1 = deg2rad(interp1(alt, long, Altitude0_1, 'pchip'));

Elev_1 = linspace(Elev_i,Elev_s,M);
Azim_1 = linspace(Azim_i,Azim_s,M);

theta_1 = pi/2 - latitude_1;
phi_1 = longitude_1;

% Radial Distance from Earth's center
R1_B = (Altitude0_1' + Re);

% Cartesian Coordinates
Rx1_I = R1_B' .* sin(theta_1) .* cos(phi_1);
Ry1_I = R1_B' .* sin(theta_1) .* sin(phi_1);
Rz1_I = R1_B' .* cos(theta_1);

R1_I = [Rx1_I'; Ry1_I'; Rz1_I'];

R1_I_mag = sqrt(Rx1_I.^2 + Ry1_I.^2 + Rz1_I.^2);
% disp(R1_I_mag);

% Earth's rotational velocity in the inertial frame
Omega_z = 2 * pi / (24 * 60 * 60);  % Sidereal Rotation Rate (rad/s)
Omega_x = 0;
Omega_y = 0;
Omega_I = [Omega_x; Omega_y; Omega_z];

% Velocity components in spherical coordinates
V1_r =  (2/160)*(D*R1_B');
V1_theta = 0;
V1_phi = Omega_z*sin(theta_1)*(Re+hi);

% Transformation from spherical to Cartesian coordinates
Vx1_I = V1_r .* sin(theta_1) .* cos(phi_1) + V1_theta.* cos(theta_1) .* cos(phi_1) - V1_phi.* sin(phi_1);
Vy1_I = V1_r .* sin(theta_1) .* sin(phi_1) + V1_theta.* cos(theta_1) .* sin(phi_1) + V1_phi.* cos(phi_1);
Vz1_I = V1_r .* cos(theta_1) - V1_theta.* sin(theta_1);

V1_I = [Vx1_I';Vy1_I'; Vz1_I'];

V1_I_mag = sqrt(Vx1_I.^2 + Vy1_I.^2 + Vz1_I.^2);
% disp(V1_I_mag);

latitude_1 = deg2rad(28);
longitude_1 = deg2rad(0);


% Attitude matrix
Q11 = - sin(latitude_1).*cos(longitude_1).*cos(Elev_1).*cos(Azim_1) - sin(longitude_1).*cos(Elev_1).*sin(Azim_1) + cos(latitude_1).*cos(longitude_1).*sin(Elev_1);

Q12 = sin(latitude_1).*cos(longitude_1).*sin(Azim_1) - sin(longitude_1).*cos(Azim_1);

Q13 = - sin(latitude_1).*cos(longitude_1).*sin(Elev_1).*cos(Azim_1) - sin(longitude_1).*sin(Elev_1).*sin(Azim_1) - cos(latitude_1).*cos(longitude_1).*cos(Elev_1);

Q21 = - sin(latitude_1).*sin(longitude_1).*cos(Elev_1).*cos(Azim_1) + cos(longitude_1).*cos(Elev_1).*sin(Azim_1) + cos(latitude_1).*sin(longitude_1).*sin(Elev_1);

Q22 = sin(latitude_1).*sin(longitude_1).*sin(Azim_1) + cos(longitude_1).*cos(Azim_1);

Q23 = - sin(latitude_1).*sin(longitude_1).*sin(Elev_1).*cos(Azim_1) + cos(longitude_1).*sin(Elev_1).*sin(Azim_1) - cos(latitude_1).*sin(longitude_1).*cos(Elev_1);

Q31 = cos(latitude_1).*cos(Elev_1).*cos(Azim_1) + sin(latitude_1).*sin(Elev_1);

Q32 = - cos(latitude_1).*sin(Azim_1);

Q33 = cos(latitude_1).*sin(Elev_1).*cos(Azim_1) - sin(latitude_1).*cos(Elev_1);

% Quaternion elements
q1 = -(0.5*(1 + Q11 - Q22 - Q33).^0.5);
q2 = ((Q12 + Q21)./(4*q1));
q3 = ((Q13 + Q31)./(4*q1));
q4 = ((Q23 - Q32)./(4*q1));

qn1 = [q1;q2;q3;q4];
test = sqrt(q1.^2+q2.^2+q3.^2+q4.^2);
end






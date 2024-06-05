function [R2_I,V2_I,qn2] = lat_long_elev_azi_vec2(M,problem) 

Elev_f = problem.Elev_f;
Azim_f = problem.Azim_f;
hf_f = problem.hf_f;
Vf_f = problem.Vf_f;
Elev_s = problem.Elev_s;
Azim_s = problem.Azim_s;
hf_s = problem.hf_s;
Vf_s = problem.Vf_s;
Re = problem.Re;


load alt_VS_inertial_lat.csv
lat = alt_VS_inertial_lat(:,1);
long = alt_VS_inertial_lat(:,2);
alt = alt_VS_inertial_lat(:,3)*1000;
Altitude0_2 = linspace(hf_s,hf_f,M);
Velocity0_2 = linspace(Vf_s,Vf_f,M);
Elev_2 = linspace(Elev_s,Elev_f,M);
Azim_2 = linspace(Azim_s,Azim_f,M);
latitude_2 = deg2rad(interp1(alt, lat, Altitude0_2, 'spline'));
longitude_2 = deg2rad(interp1(alt, long, Altitude0_2, 'spline'));

theta_2 = pi/2 - latitude_2;
phi_2 = longitude_2;
nodes = problem.nodes;
D = problem.D;


% Earth's rotational velocity in the inertial frame
Omega_z = 2 * pi / (24 * 60 * 60);  % Sidereal Rotation Rate (rad/s)
Omega_x = 0;
Omega_y = 0;
Omega_I = [Omega_x; Omega_y; Omega_z];

R2_B = (Altitude0_2 + Re);
Rx2_I = R2_B .* sin(theta_2) .* cos(phi_2);
Ry2_I = R2_B .* sin(theta_2) .* sin(phi_2);
Rz2_I = R2_B .* cos(theta_2);

R2_I = [Rx2_I; Ry2_I; Rz2_I];
R2_I_mag = sqrt(Rx2_I.^2 + Ry2_I.^2 + Rz2_I.^2);
% display(R2_I_mag);

% Velocity components in spherical coordinates
V2_r = 0;
V2_theta = Velocity0_2;
V2_phi = 0;

Vx2_I = V2_r.*sin(theta_2).*cos(phi_2) + V2_theta.*cos(theta_2).*cos(phi_2) - V2_phi.*sin(phi_2);
Vy2_I = V2_r.*sin(theta_2).*sin(phi_2) + V2_theta.*cos(theta_2).*sin(phi_2) + V2_phi.*cos(phi_2);
Vz2_I = V2_r.*cos(theta_2)-V2_theta.*sin(theta_2);

V2_I = [Vx2_I; Vy2_I; Vz2_I];
V2_I_mag = sqrt(Vx2_I.^2 + Vy2_I.^2 + Vz2_I.^2);
% disp(V2_I_mag);

latitude_2 = deg2rad(28);
longitude_2 = deg2rad(0);

% Attitude matrix
Q11 = - sin(latitude_2).*cos(longitude_2).*cos(Elev_2).*cos(Azim_2) - sin(longitude_2).*cos(Elev_2).*sin(Azim_2) + cos(latitude_2).*cos(longitude_2).*sin(Elev_2);

Q12 = sin(latitude_2).*cos(longitude_2).*sin(Azim_2) - sin(longitude_2).*cos(Azim_2);

Q13 = - sin(latitude_2).*cos(longitude_2).*sin(Elev_2).*cos(Azim_2) - sin(longitude_2).*sin(Elev_2).*sin(Azim_2) - cos(latitude_2).*cos(longitude_2).*cos(Elev_2);

Q21 = - sin(latitude_2).*sin(longitude_2).*cos(Elev_2).*cos(Azim_2) + cos(longitude_2).*cos(Elev_2).*sin(Azim_2) + cos(latitude_2).*sin(longitude_2).*sin(Elev_2);

Q22 = sin(latitude_2).*sin(longitude_2).*sin(Azim_2) + cos(longitude_2).*cos(Azim_2);

Q23 = - sin(latitude_2).*sin(longitude_2).*sin(Elev_2).*cos(Azim_2) + cos(longitude_2).*sin(Elev_2).*sin(Azim_2) - cos(latitude_2).*sin(longitude_2).*cos(Elev_2);

Q31 = cos(latitude_2).*cos(Elev_2).*cos(Azim_2) + sin(latitude_2).*sin(Elev_2);

Q32 = - cos(latitude_2).*sin(Azim_2);

Q33 = cos(latitude_2).*sin(Elev_2).*cos(Azim_2) - sin(latitude_2).*cos(Elev_2);

% Quaternion elements
q1 = -(0.5*(1 + Q11 - Q22 - Q33).^0.5);
q2 = ((Q12 + Q21)./(4*q1));
q3 = ((Q13 + Q31)./(4*q1));
q4 = ((Q23 - Q32)./(4*q1));

qn2 = [q1;q2;q3;q4];


test = sqrt(q1.^2+q2.^2+q3.^2+q4.^2);
end






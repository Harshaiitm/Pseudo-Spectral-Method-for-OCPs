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
nodes = problem.nodes;
D = problem.D;

stage_time = 160;
final_time = 1680;

load time_alt_lat_long.csv
time = time_alt_lat_long(10:35,1);
alt = time_alt_lat_long(10:35,2)*1000;
lat = time_alt_lat_long(10:35,3);
long = time_alt_lat_long(10:35,4);

t_2 = ((final_time-stage_time)/2).*nodes+(final_time+stage_time)/2;

Altitude0_2 = interp1(time, alt, t_2, 'pchip');
latitude_2 = deg2rad(interp1(alt, lat, Altitude0_2, 'pchip'));
longitude_2 = deg2rad(interp1(alt, long, Altitude0_2, 'pchip'));

Elev_2 = linspace(Elev_s,Elev_f,M);
Azim_2 = linspace(Azim_s,Azim_f,M);

theta_2 = pi/2 - latitude_2;
phi_2 = longitude_2;

% Earth's rotational velocity in the inertial frame
Omega_z = 2 * pi / (24 * 60 * 60);  % Sidereal Rotation Rate (rad/s)
Omega_x = 0;
Omega_y = 0;
Omega_I = [Omega_x; Omega_y; Omega_z];

R2_B = (Altitude0_2' + Re);
Rx2_I = R2_B' .* sin(theta_2) .* cos(phi_2);
Ry2_I = R2_B' .* sin(theta_2) .* sin(phi_2);
Rz2_I = R2_B' .* cos(theta_2);

R2_I = [Rx2_I'; Ry2_I'; Rz2_I'];
R2_I_mag = sqrt(Rx2_I.^2 + Ry2_I.^2 + Rz2_I.^2);
% display(R2_I_mag);

% Velocity components in spherical coordinates
V2_r = (2/(1680-160))*D*R2_B';
V2_theta = 0;
V2_phi = Omega_z*sin(theta_2)*(Re+hf_f);

Vx2_I = V2_r.*sin(theta_2).*cos(phi_2) + V2_theta.*cos(theta_2).*cos(phi_2) - V2_phi.*sin(phi_2);
Vy2_I = V2_r.*sin(theta_2).*sin(phi_2) + V2_theta.*cos(theta_2).*sin(phi_2) + V2_phi.*cos(phi_2);
Vz2_I = V2_r.*cos(theta_2)-V2_theta.*sin(theta_2);

V2_I = [Vx2_I'; Vy2_I'; Vz2_I'];
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






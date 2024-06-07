function [R2_I,V2_I,qn2] = lat_long_elev_azi_vec2(M,problem) 

Re = problem.Re;
nodes = problem.nodes;
D = problem.D;

stage_time = 160;
final_time = 1680;

load time_alt_lat_long.csv
time = time_alt_lat_long(:,1);
alt = time_alt_lat_long(:,2)*1000;
lat = time_alt_lat_long(:,3);
long = time_alt_lat_long(:,4);

t_2 = ((final_time-stage_time)/2).*nodes+(final_time+stage_time)/2;

Altitude0_2 = interp1(time, alt, t_2, 'pchip');
latitude_2 = deg2rad(interp1(alt, lat, Altitude0_2, 'pchip'));
longitude_2 = deg2rad(interp1(alt, long, Altitude0_2, 'pchip'));

theta_2 = pi/2 - latitude_2;
phi_2 = longitude_2;

% Earth's rotational velocity in the inertial frame
Omega_z = 2 * pi / (24 * 60 * 60);  % Sidereal Rotation Rate (rad/s)
Omega_x = 0;
Omega_y = 0;
Omega_I = [Omega_x; Omega_y; Omega_z];

R2_B = (Altitude0_2 + Re);
Rx2_I = R2_B .* sin(theta_2) .* cos(phi_2);
Ry2_I = R2_B .* sin(theta_2) .* sin(phi_2);
Rz2_I = R2_B .* cos(theta_2);

R2_I = [Rx2_I'; Ry2_I'; Rz2_I'];
R2_I_mag = sqrt(Rx2_I.^2 + Ry2_I.^2 + Rz2_I.^2);

Vx2_I = (2/(1680-160))*(D*Rx2_I);
Vy2_I = (2/(1680-160))*(D*Ry2_I);
Vz2_I = (2/(1680-160))*(D*Rz2_I);

V2_I = [Vx2_I'; Vy2_I'; Vz2_I'];
V2_I_mag = sqrt(Vx2_I.^2 + Vy2_I.^2 + Vz2_I.^2);

% Quaternion Calculations
q1 = zeros(1, M);
q2 = zeros(1, M);
q3 = zeros(1, M);
q4 = zeros(1, M);

for i = 1:M
    V_I2 = [Vx2_I(i); Vy2_I(i); Vz2_I(i)];
    V_B2 = [V2_I_mag(i); 0; 0];
    r(i, :) = vrrotvec(V_B2, V_I2);
    qn2 = (axang2quat(r(i, :)))';

    q1(i) = qn2(1);
    q2(i) = qn2(2);
    q3(i) = qn2(3);
    q4(i) = qn2(4);
end

qn2 = [q1;q2;q3;q4];
test = sqrt(q1.^2+q2.^2+q3.^2+q4.^2);
end






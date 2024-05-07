function [R2_I,V2_I,T2_I,qn2] = lat_long_elev_azi_vec2(M,problem) 

Elev_f = problem.Elev_f;
Azim_f = problem.Azim_f;
hf_f = problem.hf_f;
Vf_f = problem.Vf_f;
Elev_s = problem.Elev_s;
Azim_s = problem.Azim_s;
hf_s = problem.hf_s;
Vf_s = problem.Vf_s;
Re = problem.Re;
Thrust_max = problem.Thrust_max;
Thrust_max_2 = problem.Thrust_max_2;
Thrust_max_3 = problem.Thrust_max_3;

load alt_VS_inertial_lat.csv
lat = alt_VS_inertial_lat(:,1);
long = alt_VS_inertial_lat(:,2);
alt = alt_VS_inertial_lat(:,3)*1000;
Altitude0_2 = linspace(hf_s,hf_f,M);
Velocity0_2 = linspace(Vf_s,Vf_f,M);
Thrust0_2 = linspace(Thrust_max_2,Thrust_max_3,M);
Elev_2 = linspace(Elev_s,Elev_f,M);
Azim_2 = linspace(Azim_s,Azim_f,M);
latitude_2 = deg2rad(interp1(alt, lat, Altitude0_2, 'spline'));
longitude_2 = deg2rad(interp1(alt, long, Altitude0_2, 'spline'));

R2_B = Altitude0_2./sin(Elev_2);
Rx2_L = -cos(Elev_2).*cos(Azim_2).*R2_B;
Ry2_L = cos(Elev_2).*sin(Azim_2).*R2_B; 
Rz2_L = sin(Elev_2).*R2_B;

% Rx2_L =  -cos(Elev_2).*cos(Azim_2).*Rx2_B + sin(Azim_2).*Ry2_B - sin(Elev_2).*cos(Azim_2).*Rz2_B;
% Ry2_L =  -cos(Elev_2).*sin(Azim_2).*Rx2_B + cos(Azim_2).*Ry2_B + sin(Elev_2).*sin(Azim_2).*Rz2_B;
% Rz2_L =  sin(Elev_2).*Rx2_B- cos(Elev_2).*Rz2_B;

R2_L_mag = sqrt(Rx2_L.^2+Ry2_L.^2+Rz2_L.^2);

Rx2_I = sin(latitude_2).*cos(longitude_2).*Rx2_L - sin(longitude_2).*Ry2_L + cos(latitude_2).*cos(longitude_2).*(Re+Rz2_L);
Ry2_I = sin(latitude_2).*sin(longitude_2).*Rx2_L + cos(longitude_2).*Ry2_L + cos(latitude_2).*sin(longitude_2).*(Re+Rz2_L);
Rz2_I = -cos(latitude_2).*Rx2_L + sin(latitude_2).*(Re+Rz2_L);


R2_I = [Rx2_I;Ry2_I;Rz2_I];

R2_I_mag = sqrt(Rx2_I.^2+Ry2_I.^2+Rz2_I.^2);

V2_B = Velocity0_2./sin(Elev_2);
Vx2_L = -cos(Elev_2).*cos(Azim_2).*V2_B;
Vy2_L = cos(Elev_2).*sin(Azim_2).*V2_B; 
Vz2_L = sin(Elev_2).*V2_B;

% Vx2_L =  -cos(Elev_2).*cos(Azim_2).*Vx2_B + sin(Azim_2).*Vy2_B - sin(Elev_2).*cos(Azim_2).*Vz2_B;
% Vy2_L =  cos(Elev_2).*sin(Azim_2).*Vx2_B + cos(Azim_2).*Vy2_B + sin(Elev_2).*sin(Azim_2).*Vz2_B;
% Vz2_L =  sin(Elev_2).*Vx2_B- cos(Elev_2).*Vz2_B;

V2_L_mag = sqrt(Vx2_L.^2+Vy2_L.^2+Vz2_L.^2);

Omega_z = 2*pi/(24*60*60);      % Sideral Rotation Rate (rad/s)
Omega_x = 0; Omega_y = 0;
Omega_I = [Omega_x;Omega_y;Omega_z];

Vx2_I = sin(latitude_2).*cos(longitude_2).*Vx2_L - sin(longitude_2).*Vy2_L + cos(latitude_2).*cos(longitude_2).*Vz2_L + (Omega_y*Rz2_I-Omega_z*Ry2_I);
Vy2_I = sin(latitude_2).*sin(longitude_2).*Vx2_L + cos(longitude_2).*Vy2_L + cos(latitude_2).*sin(longitude_2).*Vz2_L + (Omega_z*Rx2_I-Omega_x*Rz2_I);
Vz2_I = -cos(latitude_2).*Vx2_L + sin(latitude_2).*Vz2_L+(Omega_x*Ry2_I-Omega_y*Rx2_I);

V2_I = [Vx2_I;Vy2_I;Vz2_I];
V2_I_mag = sqrt(Vx2_I.^2 + Vy2_I.^2 + Vz2_I.^2);

T2_B = Thrust0_2./sin(Elev_2);
Tx2_L = -cos(Elev_2).*cos(Azim_2).*T2_B;
Ty2_L = cos(Elev_2).*sin(Azim_2).*T2_B; 
Tz2_L = sin(Elev_2).*T2_B;

% Tx2_L =  -cos(Elev_2).*cos(Azim_2).*Tx2_B + sin(Azim_2).*Ty2_B - sin(Elev_2).*cos(Azim_2).*Tz2_B;
% Ty2_L =  cos(Elev_2).*sin(Azim_2).*Tx2_B + cos(Azim_2).*Ty2_B + sin(Elev_2).*sin(Azim_2).*Tz2_B;
% Tz2_L =  sin(Elev_2).*Tx2_B- cos(Elev_2).*Tz2_B;

T2_L_mag = sqrt(Tx2_L.^2 + Ty2_L.^2 + Tz2_L.^2);

Tx2_I = sin(latitude_2).*cos(longitude_2).*Tx2_L - sin(longitude_2).*Ty2_L + cos(latitude_2).*cos(longitude_2).*Tz2_L;
Ty2_I = sin(latitude_2).*sin(longitude_2).*Tx2_L + cos(longitude_2).*Ty2_L + cos(latitude_2).*sin(longitude_2).*Tz2_L;
Tz2_I = -cos(latitude_2).*Tx2_L + sin(latitude_2).*Tz2_L;

T2_I = [Tx2_I;Ty2_I;Tz2_I];
T2_I_mag = sqrt(Tx2_I.^2 + Ty2_I.^2 + Tz2_I.^2);

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
q1 = (0.5*(1 + Q11 - Q22 - Q33).^0.5);
q2 = ((Q12 + Q21)./(4*q1));
q3 = ((Q13 + Q31)./(4*q1));
q4 = ((Q23 - Q32)./(4*q1));

qn2 = [q1;q2;q3;q4];


test = sqrt(q1.^2+q2.^2+q3.^2+q4.^2);
end






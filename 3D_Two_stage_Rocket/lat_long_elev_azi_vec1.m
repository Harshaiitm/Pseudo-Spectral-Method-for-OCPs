function [R1_I,V1_I,T1_I,qn1] = lat_long_elev_azi_vec1(M,problem) 

Elev_i = problem.Elev_i;
Azim_i = problem.Azim_i;
hi = problem.hi;
Vi = problem.Vi;
Elev_s = problem.Elev_s;
Azim_s = problem.Azim_s;
hf_s = problem.hf_s;
Vf_s = problem.Vf_s;
Re = problem.Re;
Thrust_max = problem.Thrust_max;
Thrust_max_2 = problem.Thrust_max_2;

load alt_VS_inertial_lat.csv
lat = alt_VS_inertial_lat(:,1);
long = alt_VS_inertial_lat(:,2);
alt = alt_VS_inertial_lat(:,3)*1000;
Altitude0_1 = linspace(hi,hf_s,M);
Velocity0_1 = linspace(Vi,Vf_s,M);
Thrust0_1 = linspace(Thrust_max,Thrust_max_2,M);
Elev_1 = linspace(Elev_i,Elev_s,M);
Azim_1 = linspace(Azim_i,Azim_s,M);
latitude_1 = deg2rad(interp1(alt, lat, Altitude0_1, 'spline'));
longitude_1 = deg2rad(interp1(alt, long, Altitude0_1, 'spline'));


R1_B = (Altitude0_1+Re);
Rx1_I = R1_B.*cos(latitude_1).*cos(longitude_1);
Ry1_I = R1_B.*cos(latitude_1).*sin(longitude_1); 
Rz1_I = R1_B.*sin(latitude_1);


R1_I = [Rx1_I;Ry1_I;Rz1_I];

R1_I_mag = sqrt(Rx1_I.^2+Ry1_I.^2+Rz1_I.^2);

Omega_z = 2*pi/(24*60*60);      % Sideral Rotation Rate (rad/s)
Omega_x = 0; Omega_y = 0;
Omega_I = [Omega_x;Omega_y;Omega_z];

V1_B = Velocity0_1;
Vx1_I = V1_B.*cos(latitude_1).*cos(longitude_1) +  (Omega_y*Rz1_I-Omega_z*Ry1_I);
Vy1_I = V1_B.*cos(latitude_1).*sin(longitude_1) + (Omega_z*Rx1_I-Omega_x*Rz1_I);
Vz1_I = V1_B.*sin(latitude_1) + (Omega_x*Ry1_I-Omega_y*Rx1_I);

V1_I = [Vx1_I;Vy1_I;Vz1_I];
V1_I_mag = sqrt(Vx1_I.^2 + Vy1_I.^2 + Vz1_I.^2);

T1_B = Thrust0_1;
Tx1_I = T1_B.*cos(latitude_1).*cos(longitude_1);
Ty1_I = T1_B.*cos(latitude_1).*sin(longitude_1); 
Tz1_I = T1_B.*sin(latitude_1);


T1_I = [Tx1_I;Ty1_I;Tz1_I];
T1_I_mag = sqrt(Tx1_I.^2 + Ty1_I.^2 + Tz1_I.^2);

latitude_1 = deg2rad(28);
longitude_1 = deg2rad(0);
Elev_1 = deg2rad(90);
Azim_1 = deg2rad(90);

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
q1 = (0.5*(1 + Q11 - Q22 - Q33).^0.5);
q2 = ((Q12 + Q21)./(4*q1));
q3 = ((Q13 + Q31)./(4*q1));
q4 = ((Q23 - Q32)./(4*q1));

qn1 = [q1;q2;q3;q4];
test = sqrt(q1.^2+q2.^2+q3.^2+q4.^2);
end






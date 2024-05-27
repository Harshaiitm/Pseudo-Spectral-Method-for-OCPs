function [R_I,V_I,T_I,qn] = lat_long_elev_azi_vec(M,problem) 

Elev_i = problem.Elev_i;
Azim_i = problem.Azim_i;
hi = problem.hi;
Vi = problem.Vi;
Elev_f = problem.Elev_f;
Azim_f = problem.Azim_f;
hf_f = problem.hf_f;
Vf_f = problem.Vf_f;
Re = problem.Re;
Thrust_max = problem.Thrust_max;
Thrust_max_3 = problem.Thrust_max_3;

load alt_VS_inertial_lat.csv
lat = alt_VS_inertial_lat(:,1);
long = alt_VS_inertial_lat(:,2);
alt = alt_VS_inertial_lat(:,3)*1000;
Altitude = linspace(hi,hf_f,M);
Velocity = linspace(Vi,Vf_f,M);
Thrust = linspace(Thrust_max,Thrust_max_3,M);
Elev = linspace(Elev_i,Elev_i,M);
Azim = linspace(Azim_i,Azim_i,M);
latitude = deg2rad(interp1(alt, lat, Altitude, 'spline'));
longitude = deg2rad(interp1(alt, long, Altitude, 'spline'));


R_B = (Altitude+Re);
Rx_I = R_B.*cos(latitude).*cos(longitude);
Ry_I = R_B.*cos(latitude).*sin(longitude); 
Rz_I = R_B.*sin(latitude);


R_I = [Rx_I;Ry_I;Rz_I];

R_I_mag = sqrt(Rx_I.^2+Ry_I.^2+Rz_I.^2);

Omega_z = 2*pi/(24*60*60);      % Sideral Rotation Rate (rad/s)
Omega_x = 0; Omega_y = 0;
Omega_I = [Omega_x;Omega_y;Omega_z];

V_B = Velocity;
Vx_I = V_B.*cos(latitude).*cos(longitude);
Vy_I = V_B.*cos(latitude).*sin(longitude) + (Omega_z*Rx_I-Omega_x*Rz_I);
Vz_I = V_B.*sin(latitude);

V_I = [Vx_I;Vy_I;Vz_I];
V_I_mag = sqrt(Vx_I.^2 + Vy_I.^2 + Vz_I.^2);

T_B = Thrust;
Tx_I = T_B.*cos(latitude).*cos(longitude);
Ty_I = T_B.*cos(latitude).*sin(longitude); 
Tz_I = T_B.*sin(latitude);


T_I = [Tx_I;Ty_I;Tz_I];
T_I_mag = sqrt(Tx_I.^2 + Ty_I.^2 + Tz_I.^2);

latitude = deg2rad(28);
longitude = deg2rad(0);


% Attitude matrix
Q11 = - sin(latitude).*cos(longitude).*cos(Elev).*cos(Azim) - sin(longitude).*cos(Elev).*sin(Azim) + cos(latitude).*cos(longitude).*sin(Elev);

Q12 = sin(latitude).*cos(longitude).*sin(Azim) - sin(longitude).*cos(Azim);

Q13 = - sin(latitude).*cos(longitude).*sin(Elev).*cos(Azim) - sin(longitude).*sin(Elev).*sin(Azim) - cos(latitude).*cos(longitude).*cos(Elev);

Q21 = - sin(latitude).*sin(longitude).*cos(Elev).*cos(Azim) + cos(longitude).*cos(Elev).*sin(Azim) + cos(latitude).*sin(longitude).*sin(Elev);

Q22 = sin(latitude).*sin(longitude).*sin(Azim) + cos(longitude).*cos(Azim);

Q23 = - sin(latitude).*sin(longitude).*sin(Elev).*cos(Azim) + cos(longitude).*sin(Elev).*sin(Azim) - cos(latitude).*sin(longitude).*cos(Elev);

Q31 = cos(latitude).*cos(Elev).*cos(Azim) + sin(latitude).*sin(Elev);

Q32 = - cos(latitude).*sin(Azim);

Q33 = cos(latitude).*sin(Elev).*cos(Azim) - sin(latitude).*cos(Elev);

% Quaternion elements
q1 = -(0.5*(1 + Q11 - Q22 - Q33).^0.5);
q2 = ((Q12 + Q21)./(4*q1));
q3 = ((Q13 + Q31)./(4*q1));
q4 = ((Q23 - Q32)./(4*q1));

qn = [q1;q2;q3;q4];
test = sqrt(q1.^2+q2.^2+q3.^2+q4.^2);
end






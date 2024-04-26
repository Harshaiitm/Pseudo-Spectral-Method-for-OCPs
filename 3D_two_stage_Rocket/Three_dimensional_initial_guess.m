function [x0] = Three_dimensional_initial_guess(M,problem)

m0 = problem.m0;
Thrust_max = problem.Thrust_max;
Thrust_max_2 = problem.Thrust_max_2;
m0_2 = problem.m0_2;
Re = problem.Re;
mu = problem.mu;
mass2_f = problem.mass2_f;

lat_i = problem.lat_i;
long_i = problem.long_i;
Rx_i = problem.Rx_i;
Ry_i = problem.Ry_i;
Rz_i = problem.Rz_i;
Vx_i = problem.Vx_i;
Vy_i = problem.Vy_i;
Vz_i = problem.Vz_i;
hi = problem.hi;
Vi = problem.Vi;
Thrust_xi = problem.Thrust_xi; 
Thrust_yi = problem.Thrust_yi;
Thrust_zi = problem.Thrust_zi;

lat_s = problem.lat_s;
long_s = problem.long_s;
Rx_s = problem.Rx_s;
Ry_s = problem.Ry_s;
Rz_s = problem.Rz_s;
Vx_s = problem.Vx_s;
Vy_s = problem.Vy_s;
Vz_s = problem.Vz_s;
hf_s = problem.hf_s;
Vf_s = problem.Vf_s;
Thrust_xs = problem.Thrust_xs; 
Thrust_ys = problem.Thrust_ys;
Thrust_zs = problem.Thrust_zs;

lat_f = problem.lat_f;
long_f = problem.long_f;
Rx_f = problem.Rx_f;
Ry_f = problem.Ry_f;
Rz_f = problem.Rz_f;
Vx_f = problem.Vx_f;
Vy_f = problem.Vy_f;
Vz_f = problem.Vz_f;
hf = problem.hf;
Vf = problem.Vf;
Thrust_xf = problem.Thrust_xf; 
Thrust_yf = problem.Thrust_yf;
Thrust_zf = problem.Thrust_zf;
hs = 50000;

load alt_VS_inertial_lat.csv
lat = (alt_VS_inertial_lat(:,1));
long = (alt_VS_inertial_lat(:,2));
alt = (alt_VS_inertial_lat(:,3));
altitude0_1 = linspace((hi)/1000,(hs)/1000,M);
latitude_1 = deg2rad(interp1(alt, lat, altitude0_1, 'spline'));
longitude_1 = deg2rad(interp1(alt, long, altitude0_1, 'spline'));


Rx0_1 = (Re + altitude0_1*1000).* cos(latitude_1).* cos(longitude_1);
Ry0_1 = (Re + altitude0_1*1000).* cos(latitude_1).* sin(longitude_1);
Rz0_1 = (Re + altitude0_1*1000).* sin(latitude_1);
R0_1 = sqrt(Rx0_1.^2+Ry0_1.^2+Rz0_1.^2);

Velocity0_1 = linspace(Vi,sqrt(mu/(Re+hs)),M);
Vx0_1 = Velocity0_1 .* cos(latitude_1) .* cos(longitude_1);
Vy0_1 = Velocity0_1 .* cos(latitude_1) .* sin(longitude_1);
Vz0_1 = Velocity0_1 .* sin(latitude_1);
V0_1 = sqrt(Vx0_1.^2+Vy0_1.^2+Vz0_1.^2);

% 
% Thrust0_x1 = Thrust_max *  sin(Elev_i) * cos(Azim_i);
% Thrust0_y1 = Thrust_max *  sin(Elev_i) * sin(Azim_i);
% Thrust0_z1 = Thrust_max *  cos(Elev_i);

altitude0_2 = linspace((hs)/1000,(400000)/1000,10);
latitude_2 = deg2rad(interp1(alt, lat, altitude0_2, 'spline'));
longitude_2 = deg2rad(interp1(alt, long, altitude0_2, 'spline'));

Rx0_2 = (Re + altitude0_2*1000).* cos(latitude_2).* cos(longitude_1);
Ry0_2 = (Re + altitude0_2*1000).* cos(latitude_2).* sin(longitude_1);
Rz0_2 = (Re + altitude0_2*1000).* sin(latitude_2);
R0_2 = sqrt(Rx0_2.^2+Ry0_2.^2+Rz0_2.^2);

Velocity0_2 = linspace(sqrt(mu/(Re+hs)),sqrt(mu/(Re+hf)),M);
Vx0_2 = Velocity0_2 .* cos(latitude_2) .* cos(longitude_2);
Vy0_2 = Velocity0_2 .* cos(latitude_2) .* sin(longitude_2);
Vz0_2 = Velocity0_2 .* sin(latitude_2);
V0_2 = sqrt(Vx0_2.^2+Vy0_2.^2+Vz0_2.^2);

% 
x0(0*M+1:1*M) = Rx0_1;                          % Rx_1
x0(1*M+1:2*M) = Ry0_1;                          % Ry_1
x0(2*M+1:3*M) = Rz0_1;                          % Rz_1
x0(3*M+1:4*M) = Vx0_1;                          % Vx_1                                           
x0(4*M+1:5*M) = Vy0_1;                          % Vy_1
x0(5*M+1:6*M) = Vz0_1;                          % Vz_1    
x0(6*M+1:7*M)  = linspace(m0,m0_2,M);                           % mass_1
x0(7*M+1:8*M) = linspace(Thrust_xi,Thrust_xs,M);                % Thrust_x1                                
x0(8*M+1:9*M) = linspace(Thrust_yi,Thrust_ys,M);                % Thrust_y1
x0(9*M+1:10*M) = linspace(Thrust_zi,Thrust_zs,M);               % Thrust_z1
x0(10*M+1:11*M) = 0.5;                                          % q11
x0(11*M+1:12*M) = 0.5;                                          % q12
x0(12*M+1:13*M) = 0.5;                                          % q13
x0(13*M+1:14*M) = 0.5;                                          % q14
x0(14*M+1:15*M) = Rx0_2;                        % Rx_2
x0(15*M+1:16*M) = Ry0_2;                        % Ry_2 
x0(16*M+1:17*M) = Rz0_2;                        % Rz_2 
x0(17*M+1:18*M) = Vx0_2;                        % Vx_2 
x0(18*M+1:19*M) = Vy0_2;                        % Vy_2
x0(19*M+1:20*M) = Vz0_2;                        % Vz_2
x0(20*M+1:21*M) = linspace(m0_2,mass2_f,M);                     % mass_2
x0(21*M+1:22*M) = linspace(Thrust_xs,Thrust_xf,M);              % Thrust_x2
x0(22*M+1:23*M) = linspace(Thrust_ys,Thrust_yf,M);              % Thrust_y2
x0(23*M+1:24*M) = linspace(Thrust_zs,Thrust_zf,M);              % Thrust_z2
x0(24*M+1:25*M) = 0.5;                                          % q21
x0(25*M+1:26*M) = 0.5;                                          % q22
x0(26*M+1:27*M) = 0.5;                                          % q23
x0(27*M+1:28*M) = 0.5;                                          % q24
x0(28*M+1) = 137;                                               % stage_time
x0(28*M+2) = 1750;                                              % final_time

end
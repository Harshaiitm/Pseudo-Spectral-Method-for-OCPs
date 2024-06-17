clear all;
close all;
clc;

Re = 6371000; 
mu = 3.986012e14;               
Gamma = 1.4;
Rg = 287;
Temp0 = 288.16;
Omega_z = 0;     
Omega_x = 0; Omega_y = 0;
hi = 10;
hf_f = 400000;
Vi = 10;
Vf_f = sqrt(mu/(Re+hf_f));
t0 = 0;
final_time = 1683;

M = 10;
addpath('../PS_methods')                    % add the PS_method file directory
N = M-1;                            % Order of the polynomial
[nodes,weights] = LGL_nodes(N);     % calculate scaled node locations and weights
D=collocD(nodes);                   % segment differentiation matrix
  
time_span = ((final_time-t0)/2).*nodes+(final_time+t0)/2;

load alt_mach_time.csv
time = alt_mach_time(:,1);
Mach = alt_mach_time(:,3);
altitude = alt_mach_time(:,2)*1000;
Altitude = interp1(time,altitude,time_span,'spline');
Mach = interp1(time,Mach,time_span,'spline');

load alt_lat_lon.csv
lat = alt_lat_lon(:,2);
long = alt_lat_lon(:,3);
alt = alt_lat_lon(:,1)*1000;
latitude = deg2rad(interp1(alt,lat,Altitude,'spline'));
longitude = deg2rad(interp1(alt,long,Altitude,'spline'));

load Temp_vs_alt.csv
altitude = Temp_vs_alt(:,2)*1000;
Temp = Temp_vs_alt(:,1);
Static_Temp = interp1(altitude,Temp,Altitude,'spline');

R_B = (Altitude+Re);
Rx = R_B.*cos(latitude).*cos(longitude);
Ry = R_B.*cos(latitude).*sin(longitude); 
Rz = R_B.*sin(latitude);

Vrel = Mach.*sqrt(Gamma*Rg*Static_Temp);

theta = pi/2 - latitude;
phi = longitude;

% Velocity components in spherical coordinates
V_r = D * R_B;
V_theta = (D * theta) .* R_B;
V_phi = (D * phi) .* sin(theta) .* R_B;

% Transformation from spherical to Cartesian coordinates
Vx = V_r .* sin(theta) .* cos(phi) + V_theta.* cos(theta) .* cos(phi) - V_phi.* sin(phi) + (Omega_y*Rz-Omega_z*Ry);
Vy = V_r .* sin(theta) .* sin(phi) + V_theta.* cos(theta) .* sin(phi) + V_phi.* cos(phi) + (Omega_z*Rx-Omega_x*Rz);
Vz = V_r .* cos(theta) - V_theta.* sin(theta) + (Omega_x*Ry-Omega_y*Rx);


ceq(1:M,1) = D*Rx - ((final_time-t0)/2)*(Vx);
ceq(M+1:2*M,1) = D*Ry - ((final_time-t0)/2)*(Vy);
ceq(2*M+1:3*M,1) = D*Rz - ((final_time-t0)/2)*(Vz);

Re = 6371000; 
mu = 3.986012e14;               
Gamma = 1.4;
Rg = 287;
Temp0 = 288.16;
Omega_z = 2*pi/(24*60*60);     
Omega_x = 0; Omega_y = 0;
hi = 10;
hf_f = 400000;
Vi = 10;
Vf_f = sqrt(mu/(Re+hf_f));
t0 = 0;
final_time = 1650;

M = 10;
addpath('../PS_methods')                    % add the PS_method file directory
N = M-1;                            % Order of the polynomial
[nodes,weights] = LGL_nodes(N);     % calculate scaled node locations and weights
D=collocD(nodes);                   % segment differentiation matrix
  

load alt_VS_inertial_lat.csv
lat = alt_VS_inertial_lat(:,1);
long = alt_VS_inertial_lat(:,2);
alt = alt_VS_inertial_lat(:,3)*1000;
Altitude = linspace(hi,hf_f,M);
Velocity = linspace(Vi,Vf_f,M);
latitude = deg2rad(interp1(alt, lat, Altitude, 'spline'));
longitude = deg2rad(interp1(alt, long, Altitude, 'spline'));

R_B = (Altitude+Re);
Rx = R_B.*cos(latitude).*cos(longitude);
Ry = R_B.*cos(latitude).*sin(longitude); 
Rz = R_B.*sin(latitude);

total_time = ((final_time-t0)/2).*nodes+(final_time+t0)/2;

load mach_vs_time.csv
time = mach_vs_time(:,1);
Mach = mach_vs_time(:,2);
Mach = interp1(time,Mach,total_time,'spline');

load Temp_vs_alt.csv
altitude = Temp_vs_alt(:,2)*1000;
Temp = Temp_vs_alt(:,1);

% Interpolate
Static_Temp = interp1(altitude,Temp,Altitude,'spline');
Vrel = Mach'.*sqrt(Gamma*Rg*Static_Temp);

Vx = Vrel.*cos(latitude).*cos(longitude)+(Omega_y*Rz-Omega_z*Ry);
Vy = Vrel.*cos(latitude).*sin(longitude)+(Omega_z*Rx-Omega_x*Rz); 
Vz = Vrel.*sin(latitude)+(Omega_x*Ry-Omega_y*Rx);

ceq(1:M,1) = D*Rx' - ((final_time-t0)/2)*(Vx)';
ceq(M+1:2*M,1) = D*Ry' - ((final_time-t0)/2)*(Vy)';
ceq(2*M+1:3*M,1) = D*Rz' - ((final_time-t0)/2)*(Vz)';



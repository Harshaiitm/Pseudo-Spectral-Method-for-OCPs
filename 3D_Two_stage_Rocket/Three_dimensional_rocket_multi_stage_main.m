% Two Stage 3Dimensional Rocket Problem
% Three_dimensional_rocket_multi_stage_main.m
clc;clear all; close all;
%==============================================================================================%
%--- options ---%
% pseudospectral method
PS_method = 'LGL';                          % either LGL or CGL
M = 20;                                      % Number of collocation points
addpath('../PS_methods')                    % add the PS_method file directory

    if  strcmp(PS_method,'LGL')
        N = M-1;                            % Order of the polynomial
        [nodes,weights] = LGL_nodes(N);     % calculate scaled node locations and weights
        D=collocD(nodes);                   % segment differentiation matrix
  
    elseif  strcmp(PS_method,'CGL')
        N = M-1;                            % Order of the polynomial
        [nodes] = CGL_nodes(N);             % calculate scaled node locations and weights
        weights = CGL_weights(nodes);
        D=collocD(nodes);                   % segment differentiation matrix  
    end    
%================================================================================================================%
% Problem data    
Re = 6371000;                   % Radius of earth in meters
h_scale = 8500;             
mu = 3.986012e14;               % Gravitational parameter "GM" in m^3/s^2
Omega_z = 2*pi/(24*60*60);      % Sideral Rotation Rate (rad/s)
Omega_x = 0; Omega_y = 0;
rho0 = 1.225;                   % air density at Sea level 
g0 = 9.80665;                   % acceleration due to gravity at sea level
Gamma = 1.4;
R = 287;
Temp0 = 288.16;

problem.Gamma = Gamma;
problem.R = R;
problem.Temp0 = Temp0;

problem.Re = Re;
problem.h_scale = h_scale;
problem.mu = mu;
problem.Omega_x = Omega_x;
problem.Omega_y = Omega_y;
problem.Omega_z = Omega_z;
problem.rho0 = rho0;
problem.g0 = g0;

% Two stage Rocket (Kistler K-1) parameters
m0_1= 248950;                   % 1st stage total mass
mp1_by_m0_1 = 0.8324;           % 1st stage propellent fraction
mp1 = mp1_by_m0_1*m0_1;         % 1st stage total propellent
 
m0_2= 134040;                   % 2nd stage total mass
mp2_by_m0_2 = 0.88217;          % 2nd stage propellent fraction
mp2 = mp2_by_m0_2*m0_2;         % 2nd stage total propellent
m0 = m0_1 + m0_2;               % Total Rocket mass
mass1_i = m0;
mass1_f = m0-(mp1);
mass2_i = m0_2;
mass2_f = m0_2 - mp2;

problem.m0_1 = m0_1;
problem.mp1 = mp1;
problem.m0_2 = m0_2;
problem.m0 = m0;
problem.mass1_i = m0;
problem.mass1_f = mass1_f;
problem.mass2_i = mass2_i;
problem.mass2_f = mass2_f;

% Aerodynamic characteristics
A_ref = 61;                     % surface area
Cbx1 = -0.6;                    % Aerodynamic coefficients same for both stages
Cbx2 = -0.6;
dCby1_by_dbeta = -4.0;          
dCby2_by_dbeta = -4.0;
dCbz1_by_dalpha = -4.0;
dCbz2_by_dalpha = -4.0;

problem.A_ref = A_ref;
problem.Cbx1 = Cbx1;
problem.Cbx2 = Cbx2;
problem.dCby1_by_dbeta = dCby1_by_dbeta;
problem.dCby2_by_dbeta = dCby2_by_dbeta;
problem.dCbz1_by_dalpha = dCbz1_by_dalpha;
problem.dCbz2_by_dalpha = dCbz2_by_dalpha;


% Rocket engine definitions
T_max_by_W = 1.2;               % Thrust to weight ratio same for both stages
Isp = 300;                      % Specific Impulse (s) 
Thrust_max = T_max_by_W*m0*g0;
Thrust_max_2 = T_max_by_W*m0_2*g0;
Thrust_max_3 = 10;
Thrust_s = Thrust_max-Thrust_max_2;

problem.Isp = Isp;
problem.Thrust_max = Thrust_max;
problem.Thrust_max_2 = Thrust_max_2;
problem.Thrust_max_3 = Thrust_max_3;
problem.Thrust_s = Thrust_s;


% Trajectory constraints(loads)
q_max = 15000;                  % Dynamic pressure limit(Pa)
a_sen_max = 4;                  % Sensed acceleration limit (g's)

problem.q_max = q_max;
problem.a_sen_max = a_sen_max;

% Initial State
t0 = 0;
lat_i = deg2rad(28);
long_i = deg2rad(0);
hi = 10;
Vi = 10;
Elev_i = deg2rad(90);
Azim_i = deg2rad(90);

% Stage State
lat_s = deg2rad(27);
long_s = deg2rad(1);
hf_s = 50000;
Vf_s = sqrt(mu/(Re+hf_s));
Elev_s = deg2rad(70);
Azim_s = deg2rad(90);

% Final State
lat_f = deg2rad(-3);
long_f = deg2rad(87);
hf_f = 400000;
Vf_f = sqrt(mu/(Re+hf_f));
Elev_f = deg2rad(40);
Azim_f = deg2rad(90);
gamma_f = deg2rad(0);
inclin_f = deg2rad(28);

problem.lat_i = lat_i;
problem.long_i = long_i;
problem.Elev_i = Elev_i;
problem.Azim_i = Azim_i;
problem.hi = hi;
problem.Vi = Vi;
problem.t0 = t0;

problem.lat_s = lat_s;
problem.long_s = long_s;
problem.Elev_s = Elev_s;
problem.Azim_s = Azim_s;
problem.hf_s = hf_s;
problem.Vf_s = Vf_s;

problem.lat_f = lat_f;
problem.long_f = long_f;
problem.Elev_f = Elev_f;
problem.Azim_f = Azim_f;
problem.hf_f = hf_f;
problem.Vf_f = Vf_f;
problem.gamma_f = gamma_f;
problem.inclin_f = inclin_f;

% Decision veriables
x = zeros(28*M+2);
Rx_1 = x(0*M+1:1*M);
Ry_1 = x(1*M+1:2*M);
Rz_1 = x(2*M+1:3*M);
Vx_1 = x(3*M+1:4*M);
Vy_1 = x(4*M+1:5*M);
Vz_1 = x(5*M+1:6*M);
mass_1 = x(6*M+1:7*M);
Thrust_x1 = x(7*M+1:8*M);
Thrust_y1 = x(8*M+1:9*M);
Thrust_z1 = x(9*M+1:10*M);
q11 = x(10*M+1:11*M);
q12 = x(11*M+1:12*M);
q13 = x(12*M+1:13*M);
q14 = x(13*M+1:14*M);
Rx_2 = x(14*M+1:15*M);
Ry_2 = x(15*M+1:16*M);
Rz_2 = x(16*M+1:17*M);
Vx_2 = x(17*M+1:18*M);
Vy_2 = x(18*M+1:19*M);
Vz_2 = x(19*M+1:20*M);
mass_2 = x(20*M+1:21*M);
Thrust_x2 = x(21*M+1:22*M);
Thrust_y2 = x(22*M+1:23*M);
Thrust_z2 = x(23*M+1:24*M);
q21 = x(24*M+1:25*M);
q22 = x(25*M+1:26*M);
q23 = x(26*M+1:27*M);
q24 = x(27*M+1:28*M);
stage_time = x(28*M+1);
final_time = x(28*M+2);


% Initial guess for decision variables
x0 = Three_dimensional_initial_guess_new_vec(M,problem);


% Lower and Upper bounds for the variables
[lb, ub] = Three_dimensional_lower_upper_bounds(M,problem);

% linear Constraints
A = [];
Aeq = [];
b = [];
beq = [];

tic;
options =  optimoptions ('fmincon','Algorithm','sqp','Display','iter','OptimalityTolerance',...
1e-10 , 'stepTolerance', 1e-6, 'ConstraintTolerance' ,1e-7, 'MaxIterations',20,'MaxFunctionEvaluations',...
200000);
   
    if strcmp(PS_method,'LGL')
       [x,fval,ef,output] = fmincon(@(x) Three_dimensional_objective_func(x,M,problem),x0,A,b,Aeq,beq,lb,ub,@(x) Three_dimensional_Nonlinear_func_LGL(x,M,D,problem),options);
    elseif strcmp(PS_method,'CGL')
       [x,fval,ef,output] = fmincon(@(x) Three_dimensional_objective_func(x,M,problem),x0,A,b,Aeq,beq,lb,ub,@(x) Three_dimensional_Nonlinear_func_CGL(x,M,D,problem),options);  
    end

%===========================================================================================================================================================================================%    
% Check if x0 is within bounds
if all(x0 >= lb) && all(x0 <= ub)
    disp('Initial guess x0 is within the bounds.');
else
    disp('Initial guess x0 is not within the bounds.');
    
    % Find indices of elements of x0 that are not within bounds
    indices_outside_bounds = find(x0 < lb | x0 > ub);
    disp('Indices of elements of x0 not within bounds:');
    disp(indices_outside_bounds);
end


% Decision Variables
Rx_1 = x(0*M+1:1*M);
Ry_1 = x(1*M+1:2*M);
Rz_1 = x(2*M+1:3*M);
Vx_1 = x(3*M+1:4*M);
Vy_1 = x(4*M+1:5*M);
Vz_1 = x(5*M+1:6*M);
mass_1 = x(6*M+1:7*M);
Thrust_x1 = x(7*M+1:8*M);
Thrust_y1 = x(8*M+1:9*M);
Thrust_z1 = x(9*M+1:10*M);
q11 = x(10*M+1:11*M);
q12 = x(11*M+1:12*M);
q13 = x(12*M+1:13*M);
q14 = x(13*M+1:14*M);
Rx_2 = x(14*M+1:15*M);
Ry_2 = x(15*M+1:16*M);
Rz_2 = x(16*M+1:17*M);
Vx_2 = x(17*M+1:18*M);
Vy_2 = x(18*M+1:19*M);
Vz_2 = x(19*M+1:20*M);
mass_2 = x(20*M+1:21*M);
Thrust_x2 = x(21*M+1:22*M);
Thrust_y2 = x(22*M+1:23*M);
Thrust_z2 = x(23*M+1:24*M);
q21 = x(24*M+1:25*M);
q22 = x(25*M+1:26*M);
q23 = x(26*M+1:27*M);
q24 = x(27*M+1:28*M);
stage_time = x(28*M+1);
final_time = x(28*M+2);

% time span
t_1= ((stage_time-t0)/2).*nodes+(stage_time+t0)/2;
z_1 = t0:1:stage_time;
t_2= ((final_time-stage_time)/2).*nodes+(final_time+stage_time)/2;
z_2 = stage_time:1:final_time;

R_1 = sqrt(Rx_1.^2 + Ry_1.^2 + Rz_1.^2);
R_2 = sqrt(Rx_2.^2 + Ry_2.^2 + Rz_2.^2);

% Latitude_and_Longitude calculation
lat_1 = rad2deg(asin(Rz_1./R_1));
long_1 = rad2deg(acos(Rx_1./(R_1.*cos(deg2rad(lat_1)))));
% Latitude_and_Longitude calculation
lat_2 = rad2deg(asin(Rz_2./R_2));
long_2 = rad2deg(acos(Rx_2./(R_2.*cos(deg2rad(lat_2)))));

h_1 = sqrt(Rx_1.^2 + Ry_1.^2 + Rz_1.^2)-Re;
h_2 = sqrt(Rx_2.^2 + Ry_2.^2 + Rz_2.^2)-Re;



% Calculate temperature based on altitude
altitude = 0:500000;
Temp = zeros(size(altitude));
Press = zeros(size(altitude));
rho = zeros(size(altitude));

% Troposhere 
Temp(altitude < 11000) = Temp0 - 0.0065 * altitude(altitude < 11000);
Temp_Tropo =  Temp(altitude < 11000);

Press(altitude < 11000) = 101.29 * (Temp_Tropo./288.16).^5.256; 
Press_Tropo = Press(altitude < 11000)*1000;

rho(altitude < 11000) = Press_Tropo./(R*Temp_Tropo);
rho_Tropo = rho(altitude < 11000);

% Stratosphere
Temp(altitude >= 11000 & altitude < 25000) = 216.66;
Temp_Strato = Temp(altitude >= 11000 & altitude < 25000);

Press(altitude >= 11000 & altitude < 25000) = 22.65 * exp(1.73-0.000157*altitude(altitude >= 11000 & altitude < 25000)); 
Press_Strato = Press(altitude >= 11000 & altitude < 25000)*1000;

rho(altitude >= 11000 & altitude < 25000) = Press_Strato./(R*Temp_Strato);
rho_Strato = rho(altitude >= 11000 & altitude < 25000);

% Thermosphere_till_hf_s
Temp(altitude >= 25000 & altitude < hf_s) = 141.79 + 0.00299 * altitude(altitude >= 25000 & altitude < hf_s);
Temp_Thermo_1 = Temp(altitude >= 25000 & altitude < hf_s);

Press(altitude >= 25000 & altitude < hf_s) = 2.488 * (Temp_Thermo_1/216.16).^(-11.388);
Press_Thermo_1 = Press(altitude >= 25000 & altitude < hf_s)*1000; 

rho(altitude >= 25000 & altitude < hf_s) = Press_Thermo_1./(R*Temp_Thermo_1);
rho_Thermo_1 = rho(altitude >= 25000 & altitude < hf_s);

% Thermosphere_above_hf_s and till hf_f
Temp(altitude >= hf_s) = 141.79 + 0.00299 * hf_s;
Temp_Thermo_2 = Temp(altitude >= hf_s);

Press(altitude >= hf_s) = 2.488 * (Temp_Thermo_2/216.16).^(-11.388);
Press_Thermo_2 = Press(altitude >= hf_s)*1000;

rho(altitude >= hf_s) = Press_Thermo_2/(R*Temp_Thermo_2);
rho_Thermo_2 = rho(altitude >= hf_s);

% Interpolate
Static_Temp_1 = interp1(altitude, Temp,h_1);
Static_Temp_2 = interp1(altitude, Temp,h_2);
Static_Press_1 = interp1(altitude, Press,h_1);
Static_Press_2 = interp1(altitude, Press,h_2);
rho_1 = interp1(altitude,rho,h_1);
rho_2 = interp1(altitude,rho,h_2);

Vrel_x1 = Vx_1 - Rz_1.*Omega_y + Ry_1.*Omega_z;
Vrel_y1 = Vy_1 - Rx_1.*Omega_z + Rz_1.*Omega_x;
Vrel_z1 = Vz_1 - Ry_1.*Omega_x + Rx_1.*Omega_y;
Vrel_x2 = Vx_2 - Rz_2.*Omega_y + Ry_2.*Omega_z;
Vrel_y2 = Vy_2 - Rx_2.*Omega_z + Rz_2.*Omega_x;
Vrel_z2 = Vz_2 - Ry_2.*Omega_x + Rx_2.*Omega_y;

Vrel_1 = sqrt(Vrel_x1.^2 + Vrel_y1.^2 + Vrel_z1.^2);
Vrel_2 = sqrt(Vrel_x2.^2 + Vrel_y2.^2 + Vrel_z2.^2);

Mach_1 = Vrel_1./sqrt(Gamma*R*Static_Temp_1);
Mach_2 = Vrel_2./sqrt(Gamma*R*Static_Temp_2);

% rho_1 = rho0.*(Static_Temp_1/Temp0).^(Gamma-1);
% rho_2 = rho0.*(Static_Temp_2/Temp0).^(Gamma-1);
% rho_1 = rho0*exp(-(R_1-Re)./h_scale);
% rho_2 = rho0*exp(-(R_2-Re)./h_scale);

q_mag1 = 0.5* rho_1.* Vrel_1.^2;
q_mag2 = 0.5* rho_2.* Vrel_2.^2;

 
% Gravity
g_x1 = (-mu*Rx_1)./(Rx_1.^2 + Ry_1.^2 + Rz_1.^2).^(3/2);
g_y1 = (-mu*Ry_1)./(Rx_1.^2 + Ry_1.^2 + Rz_1.^2).^(3/2);
g_z1 = (-mu*Rz_1)./(Rx_1.^2 + Ry_1.^2 + Rz_1.^2).^(3/2);
g_x2 = (-mu*Rx_2)./(Rx_2.^2 + Ry_2.^2 + Rz_2.^2).^(3/2);
g_y2 = (-mu*Ry_2)./(Rx_2.^2 + Ry_2.^2 + Rz_2.^2).^(3/2);
g_z2 = (-mu*Rz_2)./(Rx_2.^2 + Ry_2.^2 + Rz_2.^2).^(3/2);

% Senced acceleration calculation
a_sen_x1 = D*Vx_1' - g_x1';
a_sen_y1 = D*Vy_1' - g_y1';
a_sen_z1 = D*Vz_1' - g_z1';
a_sen_mag1 = sqrt((a_sen_x1).^2 + (a_sen_y1).^2 + (a_sen_z1).^2);

a_sen_x2 = D*Vx_2' - g_x2';
a_sen_y2 = D*Vy_2' - g_y2';
a_sen_z2 = D*Vz_2' - g_z2';
a_sen_mag2 = sqrt((a_sen_x2).^2 + (a_sen_y2).^2 + (a_sen_z2).^2);

% Lagrange interpolation for altitude
% stage_1
collocation_points = t_1';
function_value = h_1;
altitude_1 = lagrange_interpolation_n(collocation_points, function_value, z_1);
% stage_2
collocation_points = t_2';
function_value = h_2;
altitude_2 = lagrange_interpolation_n(collocation_points, function_value, z_2);
% Multi_stage
altitude = [altitude_1 altitude_2];
z = [z_1 z_2];

% figure
figure(1)
plot(z',altitude'/1000,'g-','LineWidth',2)
xlabel('time [s]')
ylabel('Altitude [km]')
hold on
load alt_VS_time.csv
ai1 = alt_VS_time(:,1);
ai2 = alt_VS_time(:,2);
plot(ai1,ai2,'r--','LineWidth',2)
legend("PS Method","NPSOL");
title("Altitude variation w.r.t time",PS_method)
set(gca, 'FontSize', 20);
hold off
grid on

% Lagrange interpolation for Vehicle speed(Mach number)
% stage_1
collocation_points = t_1';
function_value = Mach_1;
Mach_1 = lagrange_interpolation_n(collocation_points, function_value, z_1);
% stage_2
collocation_points = t_2';
function_value = Mach_2;
Mach_2 = lagrange_interpolation_n(collocation_points, function_value, z_2);
% Multi_stage
Mach = [Mach_1 Mach_2];

figure(2)
plot(z',Mach,'g-','LineWidth',2 )
xlabel('Time [s]')
ylabel('Mach number')
hold on
load mach_vs_time.csv
al1 = mach_vs_time(:,1);
al2 = mach_vs_time(:,2);
plot(al1,al2,'r--','LineWidth',2);
grid on
legend("PS Method","NPSOL");
title("Vehicle speed w.r.t time",PS_method)
set(gca, 'FontSize', 20);
hold off 


% Lagrange interpolation for Vehicle mass 
% stage_1
collocation_points = t_1';
function_value = mass_1;
mass_1 = lagrange_interpolation_n(collocation_points, function_value, z_1);
% stage_2
collocation_points = t_2';
function_value = mass_2;
mass_2 = lagrange_interpolation_n(collocation_points, function_value, z_2);
% Multi_stage
mass = [mass_1 mass_2];

figure(3)
plot(z',mass','g-','LineWidth',2)
xlabel('Time [s]')
ylabel('mass [kg]')
grid on
hold on
load mass_VS_time.csv
al1 = mass_VS_time(:,1);
al2 = mass_VS_time(:,2);
plot(al1,al2,'r--','LineWidth',2);
grid on
legend("PS Method","NPSOL");
title("Vehicle mass variation w.r.t time",PS_method)
set(gca, 'FontSize', 20);
hold off

% Lagrange interpolation for Thrust 
Thrust_1 = sqrt(Thrust_x1.^2 + Thrust_y1.^2 + Thrust_z1.^2);
Thrust_2 = sqrt(Thrust_x2.^2 + Thrust_y2.^2 + Thrust_z2.^2);

% stage_1
collocation_points = t_1';
function_value = Thrust_1;
Thrust_1 = lagrange_interpolation_n(collocation_points, function_value, z_1);
% stage_2
collocation_points = t_2';
function_value = Thrust_2;
Thrust_2 = lagrange_interpolation_n(collocation_points, function_value, z_2);
% Multi_stage
Thrust = [Thrust_1 Thrust_2];

figure(4)
plot(z',Thrust'/1000,'g-','LineWidth',2)
xlabel('Time [s]')
ylabel('Thrust [kN]')
grid on
hold on
load Thrust_vs_time.csv
al1 = Thrust_vs_time(:,1);
al2 = Thrust_vs_time(:,2);
plot(al1,al2,'r--','LineWidth',2);
grid on
legend("PS Method","NPSOL");
title("Thrust variation w.r.t time",PS_method)
set(gca, 'FontSize', 20);
hold off 

%
latitude = [lat_1 lat_2];

figure(5)
plot(latitude',[h_1 h_2]'/1000,'g-','LineWidth',2)
xlabel('Inertual latitude(deg)')
ylabel('Altitude(Km)')
grid on
hold on
load alt_VS_inertial_lat.csv
al1 = alt_VS_inertial_lat(:,1);
al2 = alt_VS_inertial_lat(:,3);
plot(al1,al2,'r--','LineWidth',2);
grid on
legend("PS Method","NPSOL");
title("latitude variation w.r.t altitude",PS_method)
set(gca, 'FontSize', 20);

longitude = [long_1 long_2];
figure(6)
plot(longitude',[h_1 h_2]'/1000,'g-','LineWidth',2)
xlabel('Inertial longitude(deg)')
ylabel('Altitude(Km)')
grid on
hold on
load alt_VS_inertial_lat.csv
al1 = alt_VS_inertial_lat(:,2);
al2 = alt_VS_inertial_lat(:,3);
plot(al1,al2,'r--','LineWidth',2);
grid on
legend("PS Method","NPSOL");
title("longitude variation w.r.t altitude",PS_method)
set(gca, 'FontSize', 20);

% stage_1
collocation_points = t_1';
function_value = q_mag1;
q_mag1 = lagrange_interpolation_n(collocation_points, function_value, z_1);
% stage_2
collocation_points = t_2';
function_value = q_mag2;
q_mag2 = lagrange_interpolation_n(collocation_points, function_value, z_2);
% Multi_stage
q_mag = [q_mag1 q_mag2];

figure(7)
plot(z,q_mag/1000,'g-',"LineWidth",2)
xlabel('Time [s]')
ylabel('Dynamic Pressure [kPa]')
grid on
hold on
load dynamic_pressure_VS_time.csv
ai1 = dynamic_pressure_VS_time(:,1);
ai2 = dynamic_pressure_VS_time(:,2);
plot(ai1,ai2,'r--','LineWidth',2)
legend("PS Method","NPSOL");
title("Dynamic Pressure w.r.t time")
% set(gca, 'FontSize', 40);
hold off 

% stage_1
collocation_points = t_1';
function_value = a_sen_mag1;
a_sen_mag1 = lagrange_interpolation_n(collocation_points, function_value, z_1);
% stage_2
collocation_points = t_2';
function_value = a_sen_mag2;
a_sen_mag2 = lagrange_interpolation_n(collocation_points, function_value, z_2);
% Multi_stage
a_sen_mag = [a_sen_mag1 a_sen_mag2];

figure(8)
plot(z,a_sen_mag/g0,'g-',"LineWidth",2)
xlabel('Time [s]')
ylabel('Sensed acceleration[gs]')
grid on
hold on
% ylim([0 3]);
load sensedacce_VS_time.csv
ai1 = sensedacce_VS_time(:,1);
ai2 = sensedacce_VS_time(:,2);
plot(ai1,ai2,'r--','LineWidth',2)
legend("PS Method","NPSOL");
title("Sensed acceleration w.r.t time")
% set(gca, 'FontSize', 40);
hold off

figure(9)
plot([h_1 h_2], [rho_1 rho_2]);

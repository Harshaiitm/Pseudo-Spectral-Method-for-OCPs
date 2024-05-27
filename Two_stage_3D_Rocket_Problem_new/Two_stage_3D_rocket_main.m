% Two Stage 3Dimensional Rocket Problem
% Three_dimensional_rocket_multi_stage_main.m
clc;clear all; close all;
%==============================================================================================%
%--- options ---%
% pseudospectral method
PS_method = 'LGL';                          % either LGL or CGL
% Total number of nodes
M1 = 10;   % Nodes for the first stage
M2 = 20;   % Nodes for the second stage
M = M1 + M2; % Total nodes
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
Rg = 287;
Temp0 = 288.16;

problem.Gamma = Gamma;
problem.R = Rg;
problem.Temp0 = Temp0;

problem.Re = Re;
problem.h_scale = h_scale;
problem.mu = mu;
problem.Omega_x = Omega_x;
problem.Omega_y = Omega_y;
problem.Omega_z = Omega_z;
problem.rho0 = rho0;
problem.g0 = g0;
problem.M1 = M1;
problem.M2 = M2;
problem.nodes = nodes;

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
Cbx = -0.6;                    % Aerodynamic coefficients same for both stages
dCby_by_dbeta = -4.0;          
dCbz_by_dalpha = -4.0;


problem.A_ref = A_ref;
problem.Cbx = Cbx;
problem.dCby_by_dbeta = dCby_by_dbeta;
problem.dCbz_by_dalpha = dCbz_by_dalpha;

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
a_sen_max = 4*g0;                  % Sensed acceleration limit (g's)

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
x = zeros(15*M+2);
Rx = x(0*M+1:1*M);
Ry = x(1*M+1:2*M);
Rz = x(2*M+1:3*M);
Vx = x(3*M+1:4*M);
Vy = x(4*M+1:5*M);
Vz = x(5*M+1:6*M);
mass = x(6*M+1:7*M);
Thrust = x(7*M+1:8*M);
uTx = x(8*M+1:9*M);
uTy = x(9*M+1:10*M);
uTz = x(10*M+1:11*M);
q1 = x(11*M+1:12*M);
q2 = x(12*M+1:13*M);
q3 = x(13*M+1:14*M);
q4 = x(14*M+1:15*M);
stage_time = x(15*M+1);
final_time = x(15*M+2);

% Initial guess for decision variables
x0 = Two_stage_3D_rocket_initial_guess_new_vec(M,problem);


% Lower and Upper bounds for the variables
[lb, ub] = Two_stage_3D_rocket_lower_upper_bounds(M,problem);

% linear Constraints
A = [];
Aeq = [];
b = [];
beq = [];

tic;
options =  optimoptions ('fmincon','Algorithm','sqp','Display','iter','OptimalityTolerance',...
1e-10 , 'stepTolerance', 1e-6, 'ConstraintTolerance' ,1e-6, 'MaxIterations',30,'MaxFunctionEvaluations',...
200000);
   
    if strcmp(PS_method,'LGL')
       [x,fval,ef,output] = fmincon(@(x) Two_stage_3D_rocket_objective_func(x,M,problem),x0,A,b,Aeq,beq,lb,ub,@(x) Two_stage_3D_rocket_Nonlinear_func_LGL(x,M,D,problem),options);
    elseif strcmp(PS_method,'CGL')
       [x,fval,ef,output] = fmincon(@(x) Two_stage_3D_rocket_objective_func(x,M,problem),x0,A,b,Aeq,beq,lb,ub,@(x) Two_stage_3D_rocket_Nonlinear_func_CGL(x,M,D,problem),options);  
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

Rx = x(0*M+1:1*M);
Ry = x(1*M+1:2*M);
Rz = x(2*M+1:3*M);
Vx = x(3*M+1:4*M);
Vy = x(4*M+1:5*M);
Vz = x(5*M+1:6*M);
mass = x(6*M+1:7*M);
Thrust = x(7*M+1:8*M);
uTx = x(8*M+1:9*M);
uTy = x(9*M+1:10*M);
uTz = x(10*M+1:11*M);
q1 = x(11*M+1:12*M);
q2 = x(12*M+1:13*M);
q3 = x(13*M+1:14*M);
q4 = x(14*M+1:15*M);
stage_time = x(15*M+1);
final_time = x(15*M+2);


R = sqrt(Rx.^2 + Ry.^2 + Rz.^2);
h = R-Re;

% Latitude_and_Longitude calculation
latitude = rad2deg(asin(Rz./R));
longitude = rad2deg(acos(Rx./(R.*cos(deg2rad(latitude)))));

% Calculate temperature based on altitude
altitude = 0:500000;
Temp = zeros(size(altitude));
Press = zeros(size(altitude));
rho = zeros(size(altitude));

% Troposhere 
Temp(altitude < 11000) = Temp0 - 0.0065 * altitude(altitude < 11000);
Temp_Tropo =  Temp(altitude < 11000);

% Stratosphere
Temp(altitude >= 11000 & altitude < 25000) = 216.66;
Temp_Strato = Temp(altitude >= 11000 & altitude < 25000);

% Thermosphere_till_hf_s
Temp(altitude >= 25000 & altitude < hf_s) = 141.79 + 0.00299 * altitude(altitude >= 25000 & altitude < hf_s);
Temp_Thermo_1 = Temp(altitude >= 25000 & altitude < hf_s);

% Thermosphere_above_hf_s and till hf_f
Temp(altitude >= hf_s) = 141.79 + 0.00299 * hf_s;
Temp_Thermo_2 = Temp(altitude >= hf_s);

% Interpolate
Static_Temp = interp1(altitude, Temp,h);

% Attitude matrix for stage_1
Q11 = q1.^2 - q2.^2 - q3.^2 + q4.^2;
Q12 = 2*(q1.*q2 + q3.*q4);
Q13 = 2*(q1.*q3 - q2.*q4);
Q21 = 2*(q1.*q2 - q3.*q4);
Q22 = -q1.^2 + q2.^2 - q3.^2 + q4.^2;
Q23 = 2*(q2.*q3 + q1.*q4);
Q31 = 2*(q1.*q3 + q2.*q4);
Q32 = 2*(q2.*q3 - q1.*q4);
Q33 = -q1.^2 - q2.^2 + q3.^2 + q4.^2;


Q = [Q11 Q12 Q13; Q21 Q22 Q23; Q31 Q32 Q33];


% Inertial Thrust vector
Thrust_x = linspace(Thrust_max,mass2_f*g0*1.2,M).*(Q11.*uTx + Q12.*uTy + Q13.*uTz);
Thrust_y = linspace(Thrust_max,mass2_f*g0*1.2,M).*(Q21.*uTx + Q22.*uTy + Q23.*uTz);
Thrust_z = linspace(Thrust_max,mass2_f*g0*1.2,M).*(Q31.*uTx + Q32.*uTy + Q33.*uTz);

% Gravity
g_x = (-mu*Rx)./(Rx.^2 + Ry.^2 + Rz.^2).^(3/2);
g_y = (-mu*Ry)./(Rx.^2 + Ry.^2 + Rz.^2).^(3/2);
g_z = (-mu*Rz)./(Rx.^2 + Ry.^2 + Rz.^2).^(3/2);

Vrel_x = Vx - Rz.*Omega_y + Ry.*Omega_z;
Vrel_y = Vy - Rx.*Omega_z + Rz.*Omega_x;
Vrel_z = Vz - Ry.*Omega_x + Rx.*Omega_y;

Vrel = sqrt(Vrel_x.^2 + Vrel_y.^2 + Vrel_z.^2);

Mach = Vrel./sqrt(Gamma*Rg*Static_Temp);

rho = rho0*exp(-(h)./h_scale);

Vbx = Q11.*Vrel_x + Q21.*Vrel_y + Q31.*Vrel_z;
Vby = Q12.*Vrel_x + Q22.*Vrel_y + Q32.*Vrel_z;
Vbz = Q13.*Vrel_x + Q23.*Vrel_y + Q33.*Vrel_z;

Vb = sqrt(Vbx.^2 + Vby.^2 + Vbz.^2);

q_mag = 0.5* rho.* Vb.^2;

% Inertial Aerodynamic Coefficients
alpha = atan(Vbz./Vbx);
beta = atan(Vby./sqrt(Vbx.^2 + Vbz.^2));
phi = atan(Vby./Vbx);

Cby = dCby_by_dbeta*beta;
Cbz = dCbz_by_dalpha*alpha;

Cx = Q11.*Cbx + Q12.*Cby + Q13.*Cbz;
Cy = Q21.*Cbx + Q22.*Cby + Q23.*Cbz;
Cz = Q31.*Cbx + Q32.*Cby + Q33.*Cbz;

A_x = q_mag.* Cx * A_ref;
A_y = q_mag.* Cy * A_ref;
A_z = q_mag.* Cz * A_ref;

% Senced acceleration calculation
a_sen_x = (Thrust_x + A_x)./mass;
a_sen_y = (Thrust_y + A_y)./mass;
a_sen_z = (Thrust_z + A_z)./mass;
a_sen_mag = sqrt((a_sen_x).^2 + (a_sen_y).^2 + (a_sen_z).^2);


% time span
t_1= ((stage_time-t0)/2).*nodes(1:M1)+(stage_time+t0)/2;
z_1 = t0:1:stage_time;
t_2= ((final_time-stage_time)/2).*nodes(M1+1:M)+(final_time+stage_time)/2;
z_2 = stage_time:1:final_time;

t = [t_1;t_2];
z = [z_1 z_2];

% Lagrange interpolation for altitude
% stage_1
collocation_points = t';
function_value = h';
altitude = lagrange_interpolation_n(collocation_points, function_value, z);


% figure
figure(1)
plot(t',h'/1000,'g-','LineWidth',2)
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
collocation_points = t';
function_value = Mach;
% Mach = lagrange_interpolation_n(collocation_points, function_value, z);

figure(2)
plot(t',Mach,'g-','LineWidth',2 )
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
collocation_points = t';
function_value = mass;
% mass = lagrange_interpolation_n(collocation_points, function_value, z);

figure(3)
plot(t',mass','g-','LineWidth',2)
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

% stage_1
collocation_points = t';
function_value = Thrust;
% Thrust = lagrange_interpolation_n(collocation_points, function_value, z);

figure(4)
plot(t',Thrust'/1000,'g-','LineWidth',2)
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

figure(5)
plot(latitude',h'/1000,'g-','LineWidth',2)
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

figure(6)
plot(longitude',h'/1000,'g-','LineWidth',2)
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
collocation_points = t';
function_value = q_mag;
% q_mag = lagrange_interpolation_n(collocation_points, function_value, z);

figure(7)
plot(t,q_mag/1000,'g-',"LineWidth",2)
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
collocation_points = t';
function_value = a_sen_mag;
% a_sen_mag = lagrange_interpolation_n(collocation_points, function_value, z);

figure(8)
plot(t,a_sen_mag/g0,'g-',"LineWidth",2)
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
subplot(2,1,1)
plot(h, rho);
title("Density atm model w.r.t time")
subplot(2,1,2)
plot(h, rho);
title("Density Exponential model w.r.t time")
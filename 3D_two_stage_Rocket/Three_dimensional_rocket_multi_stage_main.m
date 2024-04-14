% Two Stage 3Dimensional Rocket Problem
% Three_dimensional_rocket_multi_stage_main.m
clc;clear all; close all;
%==============================================================================================%
%--- options ---%
% pseudospectral method
PS_method = 'LGL';                          % either LGL or CGL
M = 5;                                      % Number of collocation points
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
dCbz2_by_dalpha = -0.4;

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

problem.Isp = Isp;
problem.Thrust_max = Thrust_max;
problem.Thrust_max_2 = Thrust_max_2;


% Trajectory constraints(loads)
q_max = 15000;                  % Dynamic pressure limit(Pa)
a_sen_max = 4;                  % Sensed acceleration limit (g's)

problem.q_max = q_max;
problem.a_sen_max = a_sen_max;

% Initial State
t0 = 0;
Lat_i = 28;
Long_i = 0;
hi = 10;
Vi = 10;
Elev_i = 90;
Azum_i = 90;

problem.hi = hi;
problem.Vi = Vi;
problem.t0 = t0;

% Final state
hf = 400000;
Vf = sqrt(mu/(Re+hf));
gamma_f = 0;
inclin_f = deg2rad(28);


problem.hf = hf;
problem.Vf = Vf;
problem.gamma_f = gamma_f;
problem.inclin_f = inclin_f;

% Decision veriables
x = zeros(24*M+2);
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
Rx_2 = x(10*M+1:11*M);
Ry_2 = x(11*M+1:12*M);
Rz_2 = x(12*M+1:13*M);
Vx_2 = x(13*M+1:14*M);
Vy_2 = x(14*M+1:15*M);
Vz_2 = x(15*M+1:16*M);
mass_2 = x(16*M+1:17*M);
Thrust_x2 = x(17*M+1:18*M);
Thrust_y2 = x(18*M+1:19*M);
Thrust_z2 = x(19*M+1:20*M);
q1 = x(20*M+1:21*M);
q2 = x(21*M+1:22*M);
q3 = x(22*M+1:23*M);
q4 = x(23*M+1:24*M);
stage_time = x(24*M+1);
final_time = x(24*M+2);


% Initial guess for decision variables
x0 = Three_dimensional_initial_guess(M,problem);


% Lower and Upper bounds for the variables
[lb, ub] = Three_dimensional_lower_upper_bounds(M,problem);

% linear Constraints
A = [];
Aeq = [];
b = [];
beq = [];

tic;
options =  optimoptions ('fmincon','Algorithm','sqp','Display','iter','OptimalityTolerance',...
1e-10 , 'stepTolerance', 1e-6, 'ConstraintTolerance' ,1e-7, 'MaxIterations', 2000,'MaxFunctionEvaluations',...
200000);
   
    if strcmp(PS_method,'LGL')
       [x,fval,ef,output] = fmincon(@(x) Three_dimensional_objective_func(x,M,m0,m0_2),x0,A,b,Aeq,beq,lb,ub,@(x) Three_dimensional_Nonlinear_func_LGL(x,M,D,problem),options);
    elseif strcmp(PS_method,'CGL')
       [x,fval,ef,output] = fmincon(@(x) Three_dimensional_objective_func(x,M,m0,m0_2),x0,A,b,Aeq,beq,lb,ub,@(x) Three_dimensional_Nonlinear_func_CGL(x,M,D,problem),options);  
    end

%===========================================================================================================================================================================================%    
% Decision Variables
Rx_1 = x(1:M);
Ry_1 = x(1*M+1:2*M);
Rz_1 = x(2*M+1:3*M);
Vx_1 = x(3*M+1:4*M);
Vy_1 = x(4*M+1:5*M);
Vz_1 = x(5*M+1:6*M);
mass_1 = x(6*M+1:7*M);
Thrust_x1 = x(7*M+1:8*M);
Thrust_y1 = x(8*M+1:9*M);
Thrust_z1 = x(9*M+1:10*M);
Rx_2 = x(10*M+1:11*M);
Ry_2 = x(11*M+1:12*M);
Rz_2 = x(12*M+1:13*M);
Vx_2 = x(13*M+1:14*M);
Vy_2 = x(14*M+1:15*M);
Vz_2 = x(15*M+1:16*M);
mass_2 = x(16*M+1:17*M);
Thrust_x2 = x(17*M+1:18*M);
Thrust_y2 = x(18*M+1:19*M);
Thrust_z2 = x(19*M+1:20*M);
q1 = x(20*M+1:21*M);
q2 = x(21*M+1:22*M);
q3 = x(22*M+1:23*M);
q4 = x(23*M+1:24*M);
stage_time = x(24*M+1);
final_time = x(24*M+2);


% time span
t_1= ((stage_time-t0)/2).*nodes+(stage_time+t0)/2;
z_1 = t0:0.1:stage_time;
t_2= ((final_time-stage_time)/2).*nodes+(final_time+stage_time)/2;
z_2 = stage_time:0.1:final_time;


R_1 = sqrt(Rx_1.^2 + Ry_1.^2 + Rz_1.^2)-Re;
R_2 = sqrt(Rx_2.^2 + Ry_2.^2 + Rz_2.^2)-Re;

% Mach number calculation
Gamma = 1.4;
R = 287;
Temp0 = 288.16;
Temp_1 = Temp0 - 0.0065*(R_1-Re);
Temp_2 = Temp0 - 0.0065*(R_2-Re);

Vrel_x1 = Vx_1 - Rz_1.*Omega_y + Ry_1.*Omega_z;
Vrel_y1 = Vy_1 - Rx_1.*Omega_z + Rz_1.*Omega_x;
Vrel_z1 = Vz_1 - Ry_1.*Omega_x + Rx_1.*Omega_y;
Vrel_x2 = Vx_2 - Rz_2.*Omega_y + Ry_2.*Omega_z;
Vrel_y2 = Vy_2 - Rx_2.*Omega_z + Rz_2.*Omega_x;
Vrel_z2 = Vz_2 - Ry_2.*Omega_x + Rx_2.*Omega_y;

Vrel_1 = sqrt(Vrel_x1.^2 + Vrel_y1.^2 + Vrel_z1.^2);
Vrel_2 = sqrt(Vrel_x2.^2 + Vrel_y2.^2 + Vrel_z2.^2);

Mach_1 = Vrel_1./sqrt(Gamma*R*Temp_1);
Mach_2 = Vrel_2./sqrt(Gamma*R*Temp_2);



% Lagrange interpolation for altitude
% stage_1
collocation_points = t_1';
function_value = R_1;
altitude_1 = lagrange_interpolation_n(collocation_points, function_value, z_1);
% stage_2
collocation_points = t_2';
function_value = R_2;
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
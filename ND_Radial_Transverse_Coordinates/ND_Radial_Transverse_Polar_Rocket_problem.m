% single stage 2Dimensional Rocket Problem
% two_dimensional_rocket_single_stage_main.m
clc;clear all; close all;
%==============================================================================================%
%--- options ---%
% pseudospectral method
PS_method = 'LGL';                           % either LGL or LG or LGR or CGL
M = 20;                                     % Number of collocation points
addpath('../PS_methods')                    % add the PS_method file directory

    if  strcmp(PS_method,'LGL')
        N = M-1;                            % Order of the polynomial
        [nodes,weights] = LGL_nodes(N);     % calculate scaled node locations and weights
        D=collocD(nodes);                   % segment differentiation matrix
    
    elseif  strcmp(PS_method,'CGL')
        N = M-1;                             % Order of the polynomial
        [nodes] = CGL_nodes(N);              % calculate scaled node locations and weights
        weights = CGL_weights(nodes);
        D=collocD(nodes);                    % segment differentiation matrix  
    end    
%================================================================================================================%
% Problem data    
Re = 6378145;
h_scale = 8500;
mu = 3.986e14;
m0 = 100000;
A_ref = 40;
CD = 0.6;
rho0 = 1.225;
mp0_by_m0 = 0.99;
mp0 = mp0_by_m0*m0; 
g0 = 9.80665;
Isp = 300;
hf = 400000;
q_max = 14000;
a_sen_max = 3*g0;
T_max_by_W = 1.5;
Thrust_max = T_max_by_W*m0*g0;
t0 = 0;
final_time = 650;
hi = 10;
Vi = 10;
theta_i = 0;
Omega_z = 2*pi/(24*60*60);
Vf = sqrt(mu/(hf+Re));


problem.Re = Re;
problem.h_scale = h_scale;
problem.rho0 = rho0;
problem.mu = mu;
problem.m0 = m0;
problem.A_ref = A_ref;
problem.CD = CD;
problem.g0 =g0;
problem.Isp = Isp;
problem.t0 = t0;
problem.hi = hi;
problem.Vi = Vi;
problem.theta_i = theta_i;
probelm.Omega_z = Omega_z;
problem.hf = hf;
problem.Vf = Vf;
problem.Thrust_max = Thrust_max;
problem.T_max_by_W = T_max_by_W;
problem.q_max = q_max;
problem.a_sen_max = a_sen_max;


% Decision veriables
x = zeros(7*M+1);
R = x(0*M+1:1*M);               
theta = x(1*M+1:2*M);
V_r = x(2*M+1:3*M);
V_theta = x(3*M+1:4*M);
mass = x(4*M+1:5*M);
Thrust_r = x(5*M+1:6*M);
Thrust_theta = x(6*M+1:7*M);
final_time = x(7*M+1);

n_length = 1/Re;
n_velocity = 1/sqrt(mu/Re);
n_time = n_length/n_velocity;
n_mass = 1/m0;
n_thrust = 1/(m0*g0);

t = ((final_time-t0)/2).*nodes+(final_time+t0)/2;

% Initial guess values for decision variables
x0(0*M+1:1*M) = linspace((Re+10)*n_length,(Re+hf)*n_length,M);
x0(1*M+1:2*M) = 0;
x0(2*M+1:3*M) = linspace(Vi*n_velocity,sqrt(mu/(Re+hf))*n_velocity,M);
x0(3*M+1:4*M) = linspace((Omega_z*(Re+hi))*n_velocity,sqrt(mu/(Re+hf))*n_velocity,M);
x0(4*M+1:5*M) = linspace(m0*n_mass,(m0-mp0)*n_mass,M);
x0(5*M+1:6*M) = linspace(Thrust_max*n_thrust*cos(deg2rad(45)),0,M);
x0(6*M+1:7*M) = linspace(Thrust_max*n_thrust*sin(deg2rad(45)),0,M);
x0(7*M+1) = 600*n_time;

% linear inequality and equality constraints
A = [];
b = [];
Aeq = [];
beq = [];

% Lower and Upper bounds for the variables
lb(0*M+1:1*M) = 0;
lb(1*M+1:2*M) = -inf;
lb(2*M+1:3*M) = -inf;
lb(3*M+1:4*M) = -inf;
lb(4*M+1:5*M) = (m0-mp0)*n_mass;
lb(5*M+1:6*M) = -Thrust_max*n_thrust;
lb(6*M+1:7*M) = -Thrust_max*n_thrust;
lb(7*M+1) = 0;


ub(0*M+1:1*M) = inf;
ub(1*M+1:2*M) = inf;
ub(2*M+1:3*M) = inf;
ub(3*M+1:4*M) = inf;
ub(4*M+1:5*M) = m0*n_mass;
ub(5*M+1:6*M) = Thrust_max*n_thrust;
ub(6*M+1:7*M) = Thrust_max*n_thrust;
ub(7*M+1) = 650*n_time;

tic;
options =  optimoptions ('fmincon','Algorithm','sqp','Display','iter','OptimalityTolerance',...
1e-10 , 'ConstraintTolerance' ,1e-5, 'MaxIterations', 2000,'MaxFunctionEvaluations',...
200000);
   
    if strcmp(PS_method,'LGL')
       [x,fval,ef,output] = fmincon(@(x) ND_2D_Radial_Transverse_Polar_Rocket_objective_func(x,M,m0,n_mass),x0,A,b,Aeq,beq,lb,ub,@(x) ND_2D_Radial_Transverse_Polar_Rocket_Nonlinear_func_LGL(x,M,D,problem),options);
    end
    
% Stop the timer and display the elapsed time
elapsedTime = toc;
disp(['Elapsed time: ' num2str(elapsedTime) ' seconds']);

R = x(0*M+1:1*M);               
theta = x(1*M+1:2*M);
V_r = x(2*M+1:3*M);
V_theta = x(3*M+1:4*M);
mass = x(4*M+1:5*M);
Thrust_r = x(5*M+1:6*M);
Thrust_theta = x(6*M+1:7*M);
final_time = x(7*M+1);

% Dimensionlization
R = R./n_length;               
V_r = V_r./n_velocity;
V_theta = V_theta./n_velocity;
mass = mass./n_mass;
Thrust_r = Thrust_r./n_thrust;
Thrust_theta = Thrust_theta./n_thrust;
final_time = final_time./n_time; 
t0 = t0/n_time;
Isp = Isp./n_time;

h =  R - Re;
rho = rho0 * exp(-(1/h_scale).*(h));
Thrust_mag = sqrt(Thrust_r.^2 + Thrust_theta.^2);
g = mu./(R).^2;
Drag_r = - 0.5*rho.* V_r.*sqrt(V_r.^2 + V_theta.^2)*A_ref *CD;
Drag_theta = - 0.5*rho.* V_theta.*sqrt(V_r.^2 + V_theta.^2)*A_ref *CD;
q = 0.5*rho.*(V_r.^2 + V_theta.^2);
a_sen_r = (Thrust_r + Drag_r)./mass;
a_sen_theta = (Thrust_theta + Drag_theta)./mass;
a_sen_mag = sqrt(a_sen_r.^2 + a_sen_theta.^2);

t = ((final_time-t0)/2).*nodes+(final_time+t0)/2;
altitude = h;
V_mag = sqrt(V_r.^2 + V_theta.^2);

z = t0:0.1:final_time;                  % at time in seconds

% Lagrange interpolation
collocation_points = t';
function_value = altitude;
altitude = lagrange_interpolation_n(collocation_points, function_value, z);

function_value = V_mag;
velocity = lagrange_interpolation_n(collocation_points, function_value, z);

function_value = mass;
mass = lagrange_interpolation_n(collocation_points, function_value, z);

function_value = Thrust_mag; 
Thrust = lagrange_interpolation_n(collocation_points, function_value, z);

function_value = rad2deg(theta);
theta = lagrange_interpolation_n(collocation_points, function_value, z);

function_value = q;
q = lagrange_interpolation_n(collocation_points, function_value, z);

function_value = a_sen_mag;
a_sen_mag = lagrange_interpolation_n(collocation_points, function_value, z);

%
% figure
close all;
figure(1)
plot(z,altitude/1000,'g-','LineWidth',2)
xlabel('time [s]')
ylabel('Altitude [km]')
hold on
load alt_vs_time.csv
ai1 = alt_vs_time(:,1);
ai2 = alt_vs_time(:,2);
plot(ai1,ai2,'r--','LineWidth',2)
legend("PS Method","NPSOL");
% title("Altitude variation w.r.t time")
% set(gca, 'FontSize', 20);
hold off
grid on

figure(2)
plot(z,velocity/1000,'g-','LineWidth',2)
xlabel('Time [s]')
ylabel('velocity [km/s]')
hold on
load inertial_velocity_vs_time.csv
ai1 = inertial_velocity_vs_time(:,1);
ai2 = inertial_velocity_vs_time(:,2);
plot(ai1,ai2,'r--','LineWidth',2)
legend("PS Method","NPSOL");
title("Velocity variation w.r.t time")
% set(gca, 'FontSize', 40);
hold off 

figure(3)
plot(z,mass/1000,'g-','LineWidth',2)
xlabel('Time [s]')
ylabel('mass [kg]')
grid on
hold on
load mass_VS_time.csv
ai1 = mass_VS_time(:,1);
ai2 = mass_VS_time(:,2);
plot(ai1,ai2,'r--','LineWidth',2)
legend("PS Method","NPSOL");
% title("Vehicle mass variation w.r.t time")
% set(gca, 'FontSize', 40);
hold off

figure(4)
plot(z,Thrust/1000,'g-','LineWidth',2)
xlabel('Time [s]')
ylabel('Thrust [kN]')
grid on
hold on
load thrust_VS_time.csv
ai1 = thrust_VS_time(:,1);
ai2 = thrust_VS_time(:,2);
plot(ai1,ai2,'r--','LineWidth',2)
legend("PS Method","NPSOL");
title("Thrust variation w.r.t time")
% set(gca, 'FontSize', 40);
hold off 

figure(6)
plot(theta,altitude/1000,'g-',"LineWidth",2)
xlabel('Downrange angle [degree]')
ylabel('altiude')
grid on
hold on
load alt_VS_downrange.csv
ai1 = alt_VS_downrange(:,1);
ai2 = alt_VS_downrange(:,2);
plot(ai1,ai2,'r--','LineWidth',2)
legend("PS Method","NPSOL");
title("Downrange angle w.r.t altitude")
% set(gca, 'FontSize', 40);
hold off 

 
figure(7)
plot(z,q/1000,'g-',"LineWidth",2)
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
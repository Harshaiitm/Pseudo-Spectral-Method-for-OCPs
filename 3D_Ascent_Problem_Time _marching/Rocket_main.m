clc;
clear;
close all;

tspan = [0, 1750]; % Adjusted time span for the simulation

% Problem data    
Re = 6371000;                   % Radius of earth in meters
h_scale = 8500;             
mu = 3.986012e14;               % Gravitational parameter "GM" in m^3/s^2
Omega_z = 2*pi/(24*60*60);      % Sideral Rotation Rate (rad/s)
rho0 = 1.225;                   % air density at Sea level 
g0 = 9.80665;                   % acceleration due to gravity at sea level
m0_1= 248950;                   % 1st stage total mass
m0_2= 134040;                   % 2nd stage total mass
m0 = m0_1 + m0_2;

problem.Re = Re;
problem.h_scale = h_scale;
problem.mu = mu;
problem.Omega_z = Omega_z;
problem.rho0 = rho0;
problem.g0 = g0;
problem.m0 = m0;

h_i = 10 + Re;                 % h
long_i = deg2rad(0);      % longitude
lat_i = deg2rad(28);      % lat
gamma_i = deg2rad(90);    % gamma
psi_i = deg2rad(90);      % Azim
v_i = 10;
m = m0;

% Guess
x0 = [h_i; long_i; lat_i; v_i; gamma_i; psi_i;m]; 

% ODE45 solver
Opt = odeset('Events', @Limit_event);

[t, x] = ode45(@Rocket_Dynamics, tspan, x0,Opt);   

% Corrected altitude calculation
alt = (x(:,1) - (Re))./ 1000; % Adjusted from 400000 to 10 m
vel = x(:,4)./ 1000;
long = rad2deg(x(:,2));
lat = rad2deg(x(:,3));
Azim = rad2deg(x(:,6));
m = x(:,7);

figure(1)
blue = [0.4940 0.1840 0.5560];
subplot(3,2,1)
plot(t, alt, 'color', blue, 'LineWidth', 1.5)
grid on;
xlabel('Time (s)')
ylabel('Altitude (km)')
grid minor
set(gca, 'LineWidth', 1, 'FontWeight', 'bold', 'Box', 'on');

subplot(3,2,2)
plot(t, vel, 'color', blue, 'LineWidth', 1.5)
xlabel('Time (s)')
ylabel('Velocity (km/s)')
grid on;
grid minor
set(gca, 'LineWidth', 1, 'FontWeight', 'bold', 'Box', 'on');

subplot(3,2,3)
plot(t, long, 'color', blue, 'LineWidth', 1.5)
xlabel('Time (s)')
ylabel('Longitude (deg)')
grid on;
grid minor
set(gca, 'LineWidth', 1, 'FontWeight', 'bold', 'Box', 'on');

subplot(3,2,4)
plot(t, lat, 'color', blue, 'LineWidth', 1.5)
xlabel('Time (s)')
ylabel('Latitude (deg)')
grid on;
grid minor
set(gca, 'LineWidth', 1, 'FontWeight', 'bold', 'Box', 'on');

subplot(3,2,5)
FPA = rad2deg(x(:,5));
plot(t, FPA, 'color', blue, 'LineWidth', 1.5)
xlabel('Time (s)')
ylabel('FPA (deg)')
grid on;
grid minor
set(gca, 'LineWidth', 1, 'FontWeight', 'bold', 'Box', 'on');

subplot(3,2,6)
plot(t, Azim, 'color', blue, 'LineWidth', 1.5)
xlabel('Time (s)')
ylabel('Heading Angle (deg)')
grid on;
grid minor
set(gca, 'LineWidth', 1, 'FontWeight', 'bold', 'Box', 'on');
%
figure(2)
plot(t,x(:,7))
xlabel("time(s)")
ylabel("mass")

% figure(2)
% plot(vel, alt)
% xlabel('Velocity (km/s)', 'Interpreter', 'latex')
% ylabel('Altitude (km)', 'Interpreter', 'latex')
% grid on;
% grid minor
% set(gca, 'LineWidth', 1, 'FontName', 'Arial', 'FontSize', 12, 'Box', 'on', 'Color', 'white');





clc;clear all; close all;
%==============================================================================================%
%--- options ---%
% pseudospectral method
PS_method = 'LGL';   % either LGL or LG or LGR
N =100;     % Order of the polynomial
addpath('C:\Users\Harshad\OneDrive\Desktop\min_fuel_climb\PS_methods') % add the PS_method file directory

[nodes,weights] = LGL_nodes(N); % calculate scaled node locations and weights
D=collocD(nodes); % segment differentiation matrix

Re = 6378145;
mu = 3.986e14;
S = 49.2386;
g0 = 9.80665;
Isp = 1600;
H = 7254.24;
rho0 = 1.225;
t0 = 0;
tf = 400;
t = ((tf-t0)/2).*nodes+(tf+t0)/2;


x = zeros(5*N+6);
h = x(1:N+1);
v = x(N+2:2*N+2);
gamma = x(2*N+3:3*N+3);
time = x(3*N+4);
mass = x(3*N+5:4*N+5);
alpha = x(4*N+6:5*N+6);


% % Boundary Conditions 
% h0 = 0; 
% hf = 19994.88; 
% v0 = 129.314; 
% vf = 295.092;
% gamma0 = 0; 
% gammaf = 0; 
% mass0 = 19050.864;


x0(1:N) = 0;             % altitude
x0(N+1) = 19994.88;      %  final altitude
x0(N+2:2*N+1) = 129.314;     % initial velocity
x0(2*N+2) = 295.092;     % final velocity
x0(2*N+3) = 0;           % initial gamma
x0(2*N+4:3*N+2) = 0;
x0(3*N+3) = 0;           % final gamma
x0(3*N+4) = 324;         % final time
x0(3*N+5:4*N+5) = 19050.864;    % mass
x0(4*N+6:5*N+6) = (linspace(-20,20,N+1))*pi/180;      % alpha
% x0(4*N+6:5*N+6) = 20*pi/180;      % alpha

% linear inequality and equality constraints
A = [];
b = [];
Aeq = [];
beq = [];


% Lower and Upper bounds for the variables
lb(1:N+1) = 0;
lb(N+2:2*N+2) = 5;
lb(2*N+3:3*N+3) = -40*pi/180 ;
lb(3*N+4) = 0;
lb(3*N+5:4*N+5) = 22;
lb(4*N+6:5*N+6) = -pi/4;

ub(1:N+1) = 21031.2;
ub(N+2:2*N+2) = 1000;
ub(2*N+3:3*N+3) = 40*pi/180 ;
ub(3*N+4) = 400;
ub(3*N+5:4*N+5) = 20410;
ub(4*N+6:5*N+6) = pi/4;



%==============================================================================================%
% solver fmincon
% Start the timer
tic;

options =  optimoptions ('fmincon','Algorithm','sqp','Display','iter','OptimalityTolerance',...
1e-10 , 'ConstraintTolerance' ,1e-5, 'MaxIterations', 20000,'MaxFunctionEvaluations',...
200000);
   
[x,fval,ef,output] = fmincon(@(x) climb_objective_func(x,N),x0,A,b,Aeq,beq,lb,ub,@(x) climb_Nonlinear_func(x,N,D,mu,Re,t0,g0,Isp,S),options);
   

% Stop the timer and display the elapsed time
elapsedTime = toc;
disp(['Elapsed time: ' num2str(elapsedTime) ' seconds']);

h = x(1:N+1);
v = x(N+2:2*N+2);
gamma = x(2*N+3:3*N+3);
time = x(3*N+4);
mass = x(3*N+5:4*N+5);
alpha = x(4*N+6:5*N+6);

%%
% t =linspace(1,1000,N+1);
% Lagrange interpolation
t = ((x(3*N+4)-t0)/2).*nodes+(x(3*N+4)+t0)/2;
z = t;  % at time in seconds
collocation_points=t';
function_value=h;
% syms z
altitude = lagrange_interpolation_n(collocation_points, function_value, z);
% disp(['altitude Equation:', char(Position_equation)]);
% altitude = subs(Position_equation,z,z_value);
% disp(['altitude =',char(altitude),'m']);

function_value=v;
velocity=lagrange_interpolation_n(collocation_points, function_value, z);
% disp(['Velocity Equation:', char(velocity_equation)]);
% velocity = subs(velocity_equation,z,z_value);
% disp(['Velocity =',char(velocity),'m/s']);

function_value=gamma;
acceleration=lagrange_interpolation_n(collocation_points, function_value, z);
% disp(['flight_path_angle =',char(acceleration_equation),'m/s^2']);
% flight_path_angle = subs(acceleration_equation,z,z_value);
% disp(['flight_path_angle =',char(flight_path_angle),'m/s^2']);



% figure
figure(1)
plot(v/100,altitude,'g-','LineWidth',1.5)
ylim([0 20000])
xlabel('Airspeed [m/s]')
ylabel('Altitude [m]')
hold on
load alt_VS_airspeed.csv
ai1 = alt_VS_airspeed(:,1);
ai2 = alt_VS_airspeed(:,2);
plot(ai1,ai2,'r--','LineWidth',1.5)
legend("PS Method","ICLOCS2");
title("Altitude variation w.r.t airspeed")
hold off
grid on

figure(2)
plot(t,altitude,'g-','LineWidth',1.5 )
xlim([0 tf])
ylim([0 20000])
xlabel('Time [s]')
ylabel('Altitude [m]')
hold on
load alt_vs_time.csv
al1 = alt_vs_time(:,1);
al2 = alt_vs_time(:,2);
plot(al1,al2,'r--','LineWidth',1.5);
grid on
legend("PS Method","ICLOCS2");
title("Altitude variation w.r.t time")
hold off 

figure(3)
plot(t,rad2deg(gamma),'g-','LineWidth',1.5)
xlim([0 tf])
xlabel('Time [s]')
ylabel('Flight Path Angle [deg]')
grid on
hold on
load gamma_vs_time.csv
al1 = gamma_vs_time(:,1);
al2 = gamma_vs_time(:,2);
plot(al1,al2,'r--','LineWidth',1.5);
grid on
legend("PS Method","ICLOCS2");
title("Flight Path angle variation w.r.t time")
hold off



figure(4)
plot(t,v/100,'g-','LineWidth',1.5)
xlabel('Time [s]')
ylabel('Velocity [m/s]')
grid on
hold on
load velocity_vs_time.csv
al1 = velocity_vs_time(:,1);
al2 = velocity_vs_time(:,2);
plot(al1,al2,'r--','LineWidth',1.5);
grid on
legend("PS Method","ICLOCS2");
title("Aircraft velocity variation w.r.t time")
hold off 
 
figure(5)
plot(t,mass,'g-','LineWidth',1.5)
xlabel('Time [s]')
ylabel('Aircraft Mass [kg]')
grid on
hold on
load mass_vs_time.csv
al1 = mass_vs_time(:,1);
al2 = mass_vs_time(:,2);
plot(al1,al2,'r--','LineWidth',1.5);
legend("PS Method","ICLOCS2");
title("Aircraft mass variation w.r.t time")
grid on
hold off 



figure(6)
plot(t,rad2deg(alpha),'g-','LineWidth',1.5 )
xlabel('Time [s]')
ylabel('Control Input (angle of attack) [deg]')
grid on
% ylim([0 0.1]);
hold on
load alpha_vs_time.csv
al1 = alpha_vs_time(:,1);
al2 = alpha_vs_time(:,2);
plot(al1,rad2deg(al2),'r--','LineWidth',1.5);
legend("PS Method","ICLOCS2");
title("Control input(alpha) variation w.r.t time")
hold off



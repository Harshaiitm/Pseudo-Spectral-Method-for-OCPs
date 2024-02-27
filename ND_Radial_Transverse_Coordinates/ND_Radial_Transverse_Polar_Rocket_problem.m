% single stage 2Dimensional Rocket Problem
% two_dimensional_rocket_single_stage_main.m
clc;clear all; close all;
%==============================================================================================%
%--- options ---%
% pseudospectral method
PS_method = 'LGL';   % either LGL or LG or LGR
N = 40;  % Order of the polynomial
addpath('PS_methods') % add the PS_method file directory

    if  strcmp(PS_method,'LGL')
        [nodes,weights] = LGL_nodes(N); % calculate scaled node locations and weights
        D=collocD(nodes); % segment differentiation matrix 
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
problem.hf = hf;
problem.Thrust_max = Thrust_max;
problem.T_max_by_W = T_max_by_W;
problem.q_max = q_max;
problem.a_sen_max = a_sen_max;


% Decision veriables
x = zeros(7*N+8);
R = x(1:N+1);               
theta = x(N+2:2*N+2);
V_r = x(2*N+3:3*N+3);
V_theta = x(3*N+4:4*N+4);
mass = x(4*N+5:5*N+5);
Thrust_r = x(5*N+6:6*N+6);
Thrust_theta = x(6*N+7:7*N+7);
final_time = x(7*N+8);

n_length = 1/Re;
n_velocity = 1/sqrt(mu/Re);
n_time = n_length/n_velocity;
n_mass = 1/m0;
n_thrust = 1/(m0*g0);


% Initial guess values for decision variables
x0(1:N+1) =linspace((Re+10)*n_length,(Re+hf)*n_length,N+1);
x0(N+2:2*N+2) = 0;
x0(2*N+3:3*N+3) = linspace(10*n_velocity,sqrt(mu/(Re+hf))*n_velocity,N+1);
x0(3*N+4:4*N+4) = linspace(10*n_velocity,sqrt(mu/(Re+hf))*n_velocity,N+1);
x0(4*N+5:5*N+5) = linspace(m0*n_mass,(m0-mp0)*n_mass,N+1);
x0(5*N+6:6*N+6) = linspace(Thrust_max*n_thrust,0,N+1);
x0(6*N+7:7*N+7) = linspace(Thrust_max*n_thrust,0,N+1);
x0(7*N+8) = 600*n_time;

% linear inequality and equality constraints
A = [];
b = [];
Aeq = [];
beq = [];

% Lower and Upper bounds for the variables
lb(1) = (Re+10)*n_length;
lb(2:N) = Re*n_length;
lb(N+1) = (Re+hf)*n_length;
lb(N+2) = 0;
lb(N+3:2*N+2) = 0;
lb(2*N+3) = 10*n_velocity;
lb(2*N+4:3*N+2) = 1*n_velocity;
lb(3*N+3) = sqrt(mu/(Re+hf))*n_velocity;
lb(3*N+4) = 10*n_velocity;
lb(3*N+5:4*N+3) = 1*n_velocity;
lb(4*N+4) = sqrt(mu/(Re+hf))*n_velocity;
lb(4*N+5) = m0*n_mass;
lb(4*N+6:5*N+5) = (m0-mp0)*n_mass;
lb(5*N+6:6*N+6) = 0;
lb(6*N+7:7*N+7) = 0;
lb(7*N+8) = 600*n_time;

ub(1) = (Re+10)*n_length;
ub(2:N) = 1.2*((Re+hf)*n_length);
ub(N+1) = (Re+hf)*n_length;
ub(N+2) = 0;
ub(N+3:2*N+2) = pi/2;
ub(2*N+3) = (10*n_velocity);
ub(2*N+4:3*N+2) = 1.2*sqrt(mu/((Re+hf)))*n_velocity;
ub(3*N+3) = sqrt(mu/((Re+hf)))*n_velocity;
ub(3*N+4) = (10*n_velocity);
ub(3*N+5:4*N+3) = 1.2*sqrt(mu/((Re+hf)))*n_velocity;
ub(4*N+4) = sqrt(mu/((Re+hf)))*n_velocity;
ub(4*N+5) = m0*n_mass;
ub(4*N+6:5*N+5) = m0*n_mass;
ub(5*N+6:6*N+6) = Thrust_max*n_thrust;
ub(6*N+7:7*N+7) = Thrust_max*n_thrust;
ub(7*N+8) = 1000*n_time;

tic;
options =  optimoptions ('fmincon','Algorithm','sqp','Display','iter','OptimalityTolerance',...
1e-10 , 'ConstraintTolerance' ,1e-5, 'MaxIterations', 2000,'MaxFunctionEvaluations',...
20000);
   
    if strcmp(PS_method,'LGL')
       [x,fval,ef,output] = fmincon(@(x) ND_Radial_Transverse_Polar_Rocket_objective_func(x,N,m0,n_mass),x0,A,b,Aeq,beq,lb,ub,@(x) ND_Radial_Transverse_Polar_Rocket_Nonlinear_func_LGL(x,N,D,problem),options);
    end
    
% Stop the timer and display the elapsed time
elapsedTime = toc;
disp(['Elapsed time: ' num2str(elapsedTime) ' seconds']);

R = x(1:N+1);
theta = x(N+2:2*N+2);
V_r = x(2*N+3:3*N+3);
V_theta = x(3*N+4:4*N+4);
mass = x(4*N+5:5*N+5);
Thrust_r = x(5*N+6:6*N+6);
Thrust_theta = x(6*N+7:7*N+7);
final_time = x(7*N+8);

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
Thrust_mag = (Thrust_r.^2 + Thrust_theta.^2).^0.5;
g = mu./(R).^2;
Drag_r = - 0.5*rho.* V_r.*(V_r.^2 + V_theta.^2).^0.5 *A_ref *CD;
Drag_theta = - 0.5*rho.* V_theta.*(V_r.^2 + V_theta.^2).^0.5 *A_ref *CD;
q = 0.5*rho.*(V_r.^2 + V_theta.^2);
a_sen_r = (Thrust_r + Drag_r)./mass;
a_sen_theta = (Thrust_theta + Drag_theta)./mass;
a_sen_mag = (a_sen_r.^2 + a_sen_theta.^2).^0.5;

t = ((final_time-t0)/2).*nodes+(final_time+t0)/2;
altitude = h;
velocity = (V_r.^2 + V_theta.^2).^0.5;


% figure
figure(1)
plot(t,altitude/1000,'g-','LineWidth',1.5)
xlabel('time [s]')
ylabel('Altitude [km]')
hold on
title("Altitude variation w.r.t time")
hold off
grid on

figure(2)
plot(t,velocity/1000,'g-','LineWidth',1.5 )
xlabel('Time [s]')
ylabel('velocity [km/s]')
hold on
title("Velocity variation w.r.t time")
hold off 

figure(3)
plot(t,mass,'g-','LineWidth',1.5)
xlabel('Time [s]')
ylabel('mass [kg]')
grid on
hold on
title("Vehicle mass variation w.r.t time")
hold off



figure(4)
plot(t,Thrust_mag/1000,'g-','LineWidth',1.5)
xlabel('Time [s]')
ylabel('Thrust [kN]')
grid on
hold on
title("Thrust variation w.r.t time")
hold off 
 
figure(5)
plot(t,q/1000,'g-',"LineWidth",1.5)
xlabel('Time [s]')
ylabel('Dynamic Pressure [kPa]')
grid on
hold on
title("Sensed acceleration w.r.t time")
hold off 

figure(6)
plot(t,a_sen_mag/g0,'g-',"LineWidth",1.5)
xlabel('Time [s]')
ylabel('Sensed acceleration[gs]')
grid on
hold on
title("Sensed acceleration w.r.t time")
hold off 

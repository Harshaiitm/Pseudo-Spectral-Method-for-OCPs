% single stage 2Dimensional Rocket Problem
% two_dimensional_rocket_single_stage_main.m
clc;clear all; close all;
%==============================================================================================%
%--- options ---%
% pseudospectral method
PS_method = 'LGL';   % either LGL or LG or LGR
N = 20;  % Order of the polynomial
addpath('../PS_methods') % add the PS_method file directory

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
theta_0 = 0;
gamma_0 = pi/2;
alpha_0 = 0;


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
problem.theta_0 = theta_0;
problem.gamma_0 = gamma_0;
problem.alpha_0 = alpha_0;

% Decision veriables
x = zeros(7*N+8);
R = x(1:N+1);               
theta = x(N+2:2*N+2);
V = x(2*N+3:3*N+3);
gamma = x(3*N+4:4*N+4);
mass = x(4*N+5:5*N+5);
Thrust = x(5*N+6:6*N+6);
alpha = x(6*N+7:7*N+7);
final_time = x(7*N+8);

n_length = 1/Re;
n_velocity = sqrt(Re/mu);
n_time = n_length/n_velocity;
n_mass = 1/m0;
n_force = n_mass*n_length/n_time^2;
n_theta = 1/theta_0;
n_gamma = 1/gamma_0;
n_alpha = 1/alpha_0;

% Non Dimensionlization
R = R*n_length;
theta = theta *n_theta;
V = V*n_velocity;
gamma = gamma*n_gamma;
mass = mass*n_mass;
Thrust = Thrust * n_force;
alpha = alpha * n_alpha;
final_time = final_time*n_time;

% Initial guess values for decision variables
x0(1:N+1) = linspace(1,1+hf*n_length,N+1);
x0(N+2:2*N+2) = 0;
x0(2*N+3:3*N+3) = linspace(1,sqrt(mu/(Re+hf))/10,N+1);
x0(3*N+4:4*N+4) = linspace(1,0,N+1);
x0(4*N+5:5*N+5) = linspace(1,1-mp0*n_mass,N+1);
x0(5*N+6:6*N+6) = linspace(Thrust_max*n_force,0,N+1);
x0(6*N+7:7*N+7) = 0;
x0(7*N+8) = 650*n_time;

% linear inequality and equality constraints
A = [];
b = [];
Aeq = [];
beq = [];

% Lower and Upper bounds for the variables
lb(1:N+1) = 1;
lb(N+2:2*N+2) = -2;
lb(2*N+3:3*N+3) = 1;
lb(3*N+4:4*N+4) = -2;
lb(4*N+5:5*N+5) = 1-mp0*n_mass;
lb(5*N+6:6*N+6) = 0;
lb(6*N+7:7*N+7) = -2;
lb(7*N+8) = 0;

ub(1:N+1) = inf;
ub(N+2:2*N+2) = 2;
ub(2*N+3:3*N+3) = 2*sqrt(mu/(Re+hf))/10;
ub(3*N+4:4*N+4) = 2;
ub(4*N+5:5*N+5) = 1;
ub(5*N+6:6*N+6) = Thrust_max*n_force;
ub(6*N+7:7*N+7) = 2;
ub(7*N+8) = 1000*n_time;

tic;
options =  optimoptions ('fmincon','Algorithm','sqp','Display','iter','OptimalityTolerance',...
1e-10 , 'stepTolerance', 1e-6, 'ConstraintTolerance' ,1e-5, 'MaxIterations', 2000,'MaxFunctionEvaluations',...
200000);
   
    if strcmp(PS_method,'LGL')
       [x,fval,ef,output] = fmincon(@(x) Normal_Tangential_Polar_Rocket_objective_func(x,N,m0,n_mass),x0,A,b,Aeq,beq,lb,ub,@(x) Normal_Tangential_Polar_Rocket_Nonlinear_func_LGL(x,N,D,problem),options);
    end
    
% Stop the timer and display the elapsed time
elapsedTime = toc;
disp(['Elapsed time: ' num2str(elapsedTime) ' seconds']);

R = x(1:N+1);               
theta = x(N+2:2*N+2);
V = x(2*N+3:3*N+3);
gamma = x(3*N+4:4*N+4);
mass = x(4*N+5:5*N+5);
Thrust = x(5*N+6:6*N+6);
alpha = x(6*N+7:7*N+7);
final_time = x(7*N+8);


h =  R/n_length - Re;
rho = rho0 * exp(-(1/h_scale).*(h));
g = mu./(R/n_length).^2;
g = g/g0;
g0 = 1;
Isp = Isp*n_time;
q = 0.5*rho.*(V/n_velocity.^2);
Drag = q.* A_ref *CD;
a_sen_v = (Thrust/n_force.* cos(alpha/n_alpha) - Drag)./(mass/n_mass);
a_sen_gamma = (Thrust/n_force.* sin(alpha/n_alpha))./(mass/n_mass);
a_sen_mag = sqrt(a_sen_v.^2 + a_sen_gamma.^2);
Drag = Drag * n_force;

% Dimensionlization
% Non Dimensionlization
R = R/n_length;
theta = theta /n_theta;
V = V/n_velocity;
gamma = gamma/n_gamma;
mass = mass/n_mass;
Thrust = Thrust / n_force;
alpha = alpha / n_alpha;
final_time = final_time/n_time;

t = ((final_time/n_time-t0)/2).*nodes+(final_time/n_time+t0)/2;
altitude = R-Re;
velocity = V;

%% figure
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
plot(t,Thrust/1000,'g-','LineWidth',1.5)
xlabel('Time [s]')
ylabel('Thrust [kN]')
grid on
hold on
title("Thrust variation w.r.t time")
hold off 

figure(5)
plot(t,gamma,'g-',"LineWidth",1.5)
xlabel('Flight path angle [degree]')
ylabel('altiude')
grid on
hold on
title("Flight path angle w.r.t time")
hold off 

figure(6)
plot(theta,altitude/1000,'g-',"LineWidth",1.5)
xlabel('Downrange angle [degree]')
ylabel('altiude')
grid on
hold on
title("Downrange angle w.r.t altitude")
hold off 

 
figure(7)
plot(t,q/1000,'g-',"LineWidth",1.5)
xlabel('Time [s]')
ylabel('Dynamic Pressure [kPa]')
grid on
hold on
title("Sensed acceleration w.r.t time")
hold off 

figure(8)
plot(t,a_sen_mag/g0,'g-',"LineWidth",1.5)
xlabel('Time [s]')
ylabel('Sensed acceleration[gs]')
grid on
hold on
title("Sensed acceleration w.r.t time")
hold off 


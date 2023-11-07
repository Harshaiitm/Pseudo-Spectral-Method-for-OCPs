    % single stage 2Dimensional Rocket Problem
% two_dimensional_rocket_single_stage_main.m
clc;clear all; close all;
%==============================================================================================%
%--- options ---%
% pseudospectral method
PS_method = 'LGL';   % either LGL or LG or LGR
N = 40;  % Order of the polynomial
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
n_thrust = 1/(m0*g0);


% Initial guess values for decision variables
x0(1:N+1) = linspace(1,1+hf*n_length,N+1);
x0(N+2:2*N+2) = 0;
x0(2*N+3:3*N+3) = linspace(10*n_velocity,sqrt(mu/(1+hf*n_velocity)),N+1);
x0(3*N+4:4*N+4) = linspace(0,pi/2,N+1);
x0(4*N+5:5*N+5) = linspace(1,1-mp0*n_mass,N+1);
x0(5*N+6:6*N+6) = linspace(Thrust_max*n_thrust,0,N+1);
x0(6*N+7:7*N+7) = 0;
x0(7*N+8) = 600*n_time;

% linear inequality and equality constraints
A = [];
b = [];
Aeq = [];
beq = [];

% Lower and Upper bounds for the variables
lb(1) = Re*n_length;
lb(2:N) = Re*n_length;
lb(N+1) = Re*n_length;
lb(N+2:2*N+2) = -pi/2;
lb(2*N+3) = 10*n_velocity;
lb(2*N+4:3*N+2) = (10*n_velocity);
lb(3*N+3) = sqrt(mu/(1+hf*n_length));
lb(3*N+4:4*N+4) = -pi/2;
lb(4*N+5:5*N+5) = 1-mp0*n_mass;
lb(5*N+6:6*N+6) = 0;
lb(6*N+7:7*N+7) = -pi/2;
lb(7*N+8) = 600*n_time;

ub(1) = hf*n_length;
ub(2:N) = 1.2*(hf*n_length);
ub(N+1) = hf*n_length;
ub(N+2:2*N+2) = pi/2;
ub(2*N+3) = sqrt(mu/(1+hf*n_length));
ub(2*N+4:3*N+2) = 1.2*sqrt(mu/(1+hf*n_length));
ub(3*N+3) = sqrt(mu/(1+hf*n_length));
ub(3*N+4:4*N+4) = pi/2;
ub(4*N+5:5*N+5) = m0*n_mass;
ub(5*N+6:6*N+6) = Thrust_max*n_thrust;
ub(6*N+7:7*N+7) = pi/2;
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


% Dimensionlization
R = R./n_length;               
V = V./n_velocity;
mass = mass./n_mass;
Thrust = Thrust./n_thrust;
final_time = final_time./n_time;
t0 = t0/n_time;
Isp = Isp./n_time;

h =  R - Re;
rho = rho0 * exp(-(h./h_scale));
g = mu./R.^2;
g0 = mu/Re;

q = 0.5*rho.*(V).^2;
Drag = q.* A_ref *CD;

a_sen_v = (Thrust.* cos(alpha) - Drag)./(mass);
a_sen_gamma = (Thrust.* sin(alpha))./(mass);
a_sen_mag = sqrt(a_sen_v.^2 + a_sen_gamma.^2);

t = ((final_time-t0)/2).*nodes+(final_time+t0)/2;
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


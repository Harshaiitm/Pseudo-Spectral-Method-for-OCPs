% single stage Goddard Rocket Problem
% single_stage_main.m
clc;clear all; close all;
%==============================================================================================%
%--- options ---%
% pseudospectral method
PS_method = 'LGL';   % either LGL or LG or LGR
N  = 30 ;     % Order of the polynomial
addpath('C:\Users\Harshad\OneDrive\Desktop\goddard_rocket_single_stage\PS_methods') % add the PS_method file directory

    if  strcmp(PS_method,'LGL')
        [nodes,weights] = LGL_nodes(N); % calculate scaled node locations and weights
        D=collocD(nodes);           % Phase-1 differentiation matrix

    end
%================================================================================================================%
% Problem data    
Re = 6378145;
h_scale = 8500;
mu = 3.986e14;
m0 = 5000;
m0_1 = 3000;
m0_2 = 2000;
A_ref = 10;
CD = 0.2;
rho0 = 1.225;
mp0_1 = 0.6*m0_1;
mp0_2 = 0.6*m0_2; 
g0 = 9.80665;
Isp = 300;
t0 = 0;
ts = 0;
tf = 0;
% t = ((tf-t0)/2).*nodes+(tf+t0)/2;

problem.Re = Re;
problem.h_scale = h_scale;
problem.rho0 = rho0;
problem.mu = mu;
problem.m0 = m0;
problem.m0_1 = m0_1;
problem.m0_2 = m0_2;
problem.A_ref = A_ref;
problem.CD = CD;
problem.g0 =g0;
problem.Isp = Isp;
problem.t0 = t0;
problem.ts = ts;
problem.tf = tf;


% Decision veriables
x = zeros(8*N+10);
h_1 = x(1:N+1);
v_1 = x(N+2:2*N+2);
mass_1 = x(2*N+3:3*N+3);
Thrust_1 = x(3*N+4:4*N+4);
h_2 = x(4*N+5:5*N+5);
v_2 = x(5*N+6:6*N+6);
mass_2 = x(6*N+7:7*N+7);
Thrust_2 = x(7*N+8:8*N+8);
stage_time = x(8*N+9);
final_time = x(8*N+10);

% Initial guess values for decision variables
x0(1:N+1) = 0;
x0(N+2:2*N+2) = 0;
x0(2*N+3:3*N+3) = m0_1+m0_2;
x0(3*N+4:4*N+4) = m0*g0*2;
x0(4*N+5:5*N+5) = 0;
x0(5*N+6:6*N+6) = 0;
x0(6*N+7:7*N+7) = m0_2;
x0(7*N+8:8*N+8) = m0_2*g0*2;
x0(8*N+9) = 0;      
x0(8*N+10) = 0;


% linear inequality and equality constraints
A = [];
b = [];
Aeq = [];
beq = [];

% Lower and Upper bounds for the variables

lb(1:N+1) = 0;
lb(N+2:2*N+2) = 0;
lb(2*N+3) = m0_1;
lb(2*N+4:3*N+3) = m0_1*(0.4)+m0_2;
lb(3*N+4:4*N+4) = 0;
lb(4*N+5:5*N+5) = 0;
lb(5*N+6:6*N+6) = 0;
lb(6*N+7) = m0_2;
lb(6*N+8:7*N+7) = m0_2*(0.4);
lb(7*N+8:8*N+8) = 0;
lb(8*N+9) = 0;      
lb(8*N+10) = 0;

ub(1) = inf;
ub(2:N+1) = inf;
ub(N+2) = inf;
ub(N+3:2*N+2) = inf;
ub(2*N+3:3*N+3) = m0_1+m0_2;
ub(3*N+4:4*N+4) = m0*g0*2;
ub(4*N+5:5*N+5) = inf;
ub(5*N+6:6*N+6) = inf;
ub(6*N+7:7*N+7) = m0_2;
ub(7*N+8:8*N+8) = m0_2*g0*2;
ub(8*N+9) = inf;      
ub(8*N+10) = inf;


tic;
options =  optimoptions ('fmincon','Algorithm','sqp','Display','iter','OptimalityTolerance',...
1e-10 , 'ConstraintTolerance' ,1e-5, 'MaxIterations', 2000,'MaxFunctionEvaluations',...
500000);
   
    if strcmp(PS_method,'LGL')
       [x,fval,ef,output] = fmincon(@(x) multi_stage_objective_func(x,N),x0,A,b,Aeq,beq,lb,ub,@(x) multi_stage_Nonlinear_func_LGL(x,N,D,problem),options);
    end
    
% Stop the timer and display the elapsed time
elapsedTime = toc;
disp(['Elapsed time: ' num2str(elapsedTime) ' seconds']);

h_1 = x(1:N+1);
v_1 = x(N+2:2*N+2);
mass_1 = x(2*N+3:3*N+3);
Thrust_1 = x(3*N+4:4*N+4);
h_2 = x(4*N+5:5*N+5);
v_2 = x(5*N+6:6*N+6);
mass_2 = x(6*N+7:7*N+7);
Thrust_2 = x(7*N+8:8*N+8);
stage_time = x(8*N+9);
final_time = x(8*N+10);

t(1:N+1) = ((stage_time-t0)/2).*nodes+(stage_time+t0)/2;
t(N+2:2*N+2) = ((final_time-stage_time)/2).*nodes+(final_time+stage_time)/2;

z = 0:1:final_time;  % at time in seconds


collocation_points=t';
altitude = [h_1';h_2'];
velocity = [v_1';v_2'];
mass = [mass_1';mass_2'];
Thrust = [Thrust_1';Thrust_2'];
t = [t(1:N+1)';t(N+2:2*N+2)'];
% function_value=h;
% altitude = lagrange_interpolation_n(collocation_points, function_value, z);

% function_value=v;
% velocity=lagrange_interpolation_n(collocation_points, function_value, z);

% function_value=mass;
% mass=lagrange_interpolation_n(collocation_points, function_value, z);

% function_value = Thrust;
% Thrust = lagrange_interpolation_n(collocation_points,function_value,z);

% figure
figure(1)
plot(t',altitude'/1000,'g-','LineWidth',1.5)
xlabel('time [s]')
ylabel('Altitude [km]')
hold on
load alt_VS_time.csv
ai1 = alt_VS_time(:,1);
ai2 = alt_VS_time(:,2);
plot(ai1,ai2,'r--','LineWidth',1.5)
legend("PS Method","MIT");
title("Altitude variation w.r.t time")
hold off
grid on

figure(2)
plot(t',velocity'/1000,'g-','LineWidth',1.5 )
xlabel('Time [s]')
ylabel('velocity [km/s]')
hold on
load velocity_vs_time.csv
al1 = velocity_vs_time(:,1);
al2 = velocity_vs_time(:,2);
plot(al1,al2,'r--','LineWidth',1.5);
grid on
legend("PS Method","MIT");
title("Velocity variation w.r.t time")
hold off 

figure(3)
plot(t',mass','g-','LineWidth',1.5)
xlabel('Time [s]')
ylabel('mass [kg]')
grid on
hold on
load mass_vs_time.csv
al1 = mass_vs_time(:,1);
al2 = mass_vs_time(:,2);
plot(al1,al2,'r--','LineWidth',1.5);
grid on
legend("PS Method","MIT");
title("Vehicle mass variation w.r.t time")
hold off



figure(4)
plot(t',Thrust'/1000,'g-','LineWidth',1.5)
xlabel('Time [s]')
ylabel('Thrust [kN]')
grid on
hold on
load Thrust_vs_time.csv
al1 = Thrust_vs_time(:,1);
al2 = Thrust_vs_time(:,2);
plot(al1,al2,'r--','LineWidth',1.5);
grid on
legend("PS Method","ICLOCS2");
title("Thrust variation w.r.t time")
hold off 
 


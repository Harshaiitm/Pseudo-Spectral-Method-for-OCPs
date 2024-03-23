% single stage 2Dimensional Rocket Problem
% two_dimensional_rocket_single_stage_main.m
clc;clear all; close all;
%==============================================================================================%
%--- options ---%
% pseudospectral method
PS_method = 'LGL';                           % either LGL or LG or LGR or CGL
M = 30;                                     % Number of collocation points
addpath('../PS_methods')                    % add the PS_method file directory

    if  strcmp(PS_method,'LGL')
        N = M-1;                            % Order of the polynomial
        [nodes,weights] = LGL_nodes(N);     % calculate scaled node locations and weights
        D=collocD(nodes);                   % segment differentiation matrix
    
    elseif strcmp(PS_method,'LGR')
        N = M-1;                            % Order of the polynomial
        [nodes,weights] = LGR_nodes(N);     % LGR_nodes gives the N+1 nodes in [-1 1)
        nodes = flip(-nodes);               % Flipped LGR method
        weights = flip(weights);            % weights are flipped
        nodes = [-1;nodes];                 % Introducing non-collocated point -1 
        D = collocD(nodes);                 % differentiation matrix of size M by M
        D(1,:) = [];                        % deletion of first row associated with non-collocated point
        nodes(1) = [];
    
    elseif strcmp(PS_method,'LG')
        N = M;                              % Order of the polynomial
        [nodes,weights]=LG_nodes(N,-1,1);   % calculate scaled node locations and weights
        nodes = [-1;nodes;1];                 % Introducing non-collocated point -1
        D=collocD(nodes);                   % segment differentiation matrix
        D(1,:) = [];                        % deletion of first row associated with non-collocated point
        D(end,:) = [];
        nodes(1) = [];
        nodes(end) = [];
    
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

% Decision veriables
x = zeros(7*M+1);
R = x(1:M);               
theta = x(M+1:2*M);
V = x(2*M+1:3*M);
gamma = x(3*M+1:4*M);
mass = x(4*M+1:5*M);
Thrust = x(5*M+1:6*M);
alpha = x(6*M+1:7*M);
final_time = x(7*M+1);

% Initial guess values for decision variables
x0(1:M) = linspace(Re,Re+hf,M);
x0(M+1:2*M) = 0;
x0(2*M+1:3*M) = linspace(10,sqrt(mu/(Re+hf)),M);
x0(3*M+1:4*M) = linspace(pi/2,0,M);
x0(4*M+1:5*M) = linspace(m0,m0-mp0,M);
x0(5*M+1:6*M) = linspace(Thrust_max,0,M);
x0(6*M+1:7*M) = 0;
x0(7*M+1) = 650;


% Initial and Final conditions
Ri = 1;
Vi = 1;
Vf = 1;
theta_i = 0;
gamma_i = pi/2;
gamma_f = 0;
alpha_i = 0;
mass_i = 1;



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
problem.Ri = Ri;
problem.Vi = Vi; 
problem.Vf = Vf; 
problem.theta_i = theta_i; 
problem.gamma_i = gamma_i; 
problem.gamma_f = gamma_f; 
problem.alpha_i = alpha_i; 
problem.mass_i = mass_i;


% linear inequality and equality constraints
A = [];
b = [];
Aeq = [];
beq = [];

% Lower and Upper bounds for the variables
lb(1:M) = Re;
lb(M+1:2*M) = -pi;
lb(2*M+1:3*M) = 10;
lb(3*M+1:4*M) = -pi;
lb(4*M+1:5*M) = m0-mp0;
lb(5*M+1:6*M) = 0;
lb(6*M+1:7*M) = -pi;
lb(7*M+1) = 0;

ub(1:M) = 2*(Re+hf);
ub(M+1:2*M) = pi;
ub(2*M+1:3*M) = 2*sqrt(mu/(Re+hf));
ub(3*M+1:4*M) = pi;
ub(4*M+1:5*M) = m0;
ub(5*M+1:6*M) = Thrust_max;
ub(6*M+1:7*M) = pi;
ub(7*M+1) = 1000;

tic;
options =  optimoptions ('fmincon','Algorithm','sqp','Display','iter','OptimalityTolerance',...
1e-10 , 'stepTolerance', 1e-10, 'ConstraintTolerance' ,1e-5, 'MaxIterations', 2000,'MaxFunctionEvaluations',...
20000);
   
    if strcmp(PS_method,'LGL')
       [x,fval,ef,output] = fmincon(@(x) Normal_Tangential_Polar_Rocket_objective_func(x,M,m0),x0,A,b,Aeq,beq,lb,ub,@(x) Normal_Tangential_Polar_Rocket_Nonlinear_func_LGL(x,M,D,problem),options);
    elseif strcmp(PS_method,'LGR')
       [x,fval,ef,output] = fmincon(@(x) Normal_Tangential_Polar_Rocket_objective_func(x,M,m0),x0,A,b,Aeq,beq,lb,ub,@(x) Normal_Tangential_Polar_Rocket_Nonlinear_func_LGR(x,M,D,problem),options);
    elseif strcmp(PS_method,'LG')
       [x,fval,ef,output] = fmincon(@(x) Normal_Tangential_Polar_Rocket_objective_func(x,M,m0),x0,A,b,Aeq,beq,lb,ub,@(x) Normal_Tangential_Polar_Rocket_Nonlinear_func_LG(x,M,D,problem),options);
    elseif strcmp(PS_method,'CGL')
       [x,fval,ef,output] = fmincon(@(x) Normal_Tangential_Polar_Rocket_objective_func(x,M,m0),x0,A,b,Aeq,beq,lb,ub,@(x) Normal_Tangential_Polar_Rocket_Nonlinear_func_CGL(x,M,D,problem),options);
    end
    
% Stop the timer and display the elapsed time
elapsedTime = toc;
disp(['Elapsed time: ' num2str(elapsedTime) ' seconds']);

R = x(1:M);               
theta = x(M+1:2*M);
V = x(2*M+1:3*M);
gamma = x(3*M+1:4*M);
mass = x(4*M+1:5*M);
Thrust = x(5*M+1:6*M);
alpha = x(6*M+1:7*M);
final_time = x(7*M+1);


h =  R - Re;
rho = rho0 * exp(-(1/h_scale).*(h));
g = mu./(R).^2;
q = 0.5*rho.*(V.^2);
Drag = q.* A_ref *CD;
a_sen_v = (Thrust.* cos(alpha) - Drag)./mass;
a_sen_gamma = (Thrust.* sin(alpha))./mass;
a_sen_mag = (a_sen_v.^2 + a_sen_gamma.^2).^0.5;

t = ((final_time-t0)/2).*nodes+(final_time+t0)/2;
altitude = R-Re;
velocity = V;

% Lagrange interpolation
% z = 0:0.1:x(4*N+5);  % at time in seconds
% 
% 
% collocation_points = t';
% function_value = R;
% altitude = lagrange_interpolation_n(collocation_points, function_value, z);
% 
% function_value = V;
% velocity=lagrange_interpolation_n(collocation_points, function_value, z);
% 
% function_value = mass;
% mass=lagrange_interpolation_n(collocation_points, function_value, z);
% 
% function_value = Thrust;
% Thrust = lagrange_interpolation_n(collocation_points,function_value,z);
%%
% figure
figure(1)
plot(t,altitude/1000,'g-','LineWidth',1.5)
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
plot(t,velocity/1000,'g-','LineWidth',1.5 )
xlabel('Time [s]')
ylabel('Inertial_Velocity [km/s]')
hold on
load inertial_velocity_vs_time.csv
al1 = inertial_velocity_vs_time(:,1);
al2 = inertial_velocity_vs_time(:,2);
plot(al1,al2,'r--','LineWidth',1.5);
grid on
legend("PS Method","MIT");
title("Inertial velocity variation w.r.t time")
hold off 

figure(3)
plot(t,mass/1000,'g-','LineWidth',1.5)
xlabel('Time [s]')
ylabel('mass [kg]')
grid on
hold on
load mass_vs_time.csv
al1 = mass_vs_time(:,1);
al2 = mass_vs_time(:,2);
plot(al1,al2,'r--','LineWidth',1.5);
grid on
xlim([0,final_time]);
legend("PS Method","MIT");
title("Vehicle mass variation w.r.t time")
hold off

figure(4)
plot(t,Thrust/1000,'g-','LineWidth',1.5)
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

figure(5)
plot(t,gamma,'g-',"LineWidth",1.5)
xlabel('Flight path angle [degree]')
ylabel('altiude')
grid on
hold on
load flightpath_VS_time.csv
al1 = flightpath_VS_time(:,1);
al2 = flightpath_VS_time(:,2);
plot(al1,al2,'r--','LineWidth',1.5);
grid on
legend("PS Method","ICLOCS2");
title("Flight path angle w.r.t time")
hold off 

figure(7)
plot(theta,altitude/1000,'g-',"LineWidth",1.5)
xlabel('Downrange angle [degree]')
ylabel('altiude')
grid on
hold on
load alt_VS_downrange.csv
al1 = alt_VS_downrange(:,1);
al2 = alt_VS_downrange(:,2);
plot(al1,al2,'r--','LineWidth',1.5);
grid on
legend("PS Method","ICLOCS2");
title("Downrange angle w.r.t altitude")
hold off 

 
figure(5)
plot(t,q/1000,'g-',"LineWidth",1.5)
xlabel('Time [s]')
ylabel('Dynamic Pressure [kPa]')
grid on
hold on
load dynamic_pressure_VS_time.csv
al1 = dynamic_pressure_VS_time(:,1);
al2 = dynamic_pressure_VS_time(:,2);
plot(al1,al2,'r--','LineWidth',1.5);
grid on
legend("PS Method","ICLOCS2");
title("Sensed acceleration w.r.t time")
hold off 

figure(6)
plot(t,a_sen_mag/g0,'g-',"LineWidth",1.5)
xlabel('Time [s]')
ylabel('Sensed acceleration[gs]')
grid on
hold on
load sensedacce_VS_time.csv
al1 = sensedacce_VS_time(:,1);
al2 = sensedacce_VS_time(:,2);
plot(al1,al2,'r--','LineWidth',1.5);
grid on
legend("PS Method","ICLOCS2");
title("Sensed acceleration w.r.t time")
hold off 


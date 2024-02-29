% single stage Goddard Rocket Problem
% single_stage_main.m
clc;clear all; close all;
%==============================================================================================%
%--- options ---%
% pseudospectral method
PS_method = 'LGL';   % either LGL or LG or LGR
N = 50;     % Order of the polynomial
addpath('../PS_methods') % add the PS_method file directory

    if  strcmp(PS_method,'LGL')
        [nodes,weights] = LGL_nodes(N); % calculate scaled node locations and weights
        D=collocD(nodes); % segment differentiation matrix
    elseif strcmp(PS_method,'LG')
        nodes(1) =-1;
        [nodes(2:N+1,1),weights]=LG_nodes(N,-1,1); % calculate scaled node locations and weights
        nodes(2:N+1,1) = flip(nodes(2:N+1,1)); 
        D=collocD(nodes); % segment differentiation matrix
        D(N+1,:) = [];
        nodes(N+2) = 1;
    elseif strcmp(PS_method,'LGR')
        [nodes,weights] = LGR_nodes(N); % calculate scaled node locations and weights
        nodes = flip(-nodes);
        nodes(1) = -1;
        D = collocD(nodes); % segment differentiation matrix
        weights = flip(weights);
        D(N+1,:) = [];    
    else strcmp(PS_method,'CGL')
         [nodes] = CGL_nodes(N);     % calculate scaled node locations and weights
         weights = CGL_weights(nodes);
         D=collocD(nodes);           % segment differentiation matrix  
    end   
%================================================================================================================%
% Problem data    
Re = 6378145;
h_scale = 8500;
mu = 3.986e14;
m0 = 5000;
A_ref = 10;
CD = 0.2;
rho0 = 1.225;
mp0 = 0.6*m0; 
g0 = 9.80665;
Isp = 300;
% Thrust = m0 * g0 *2;
t0 = 0;
tf = 0;
t = ((tf-t0)/2).*nodes+(tf+t0)/2;

problem.Re = Re;
problem.h_scale = h_scale;
problem.rho0 = rho0;
problem.mu = mu;
problem.m0 = m0;
problem.A_ref = A_ref;
problem.CD = CD;
problem.g0 =g0;
problem.Isp = Isp;
% problem.Thrust = Thrust;
problem.t0 = t0;
problem.tf = tf;


% Decision veriables
x = zeros(4*N+5);
h = x(1:N+1);
v = x(N+2:2*N+2);
mass = x(2*N+3:3*N+3);
Thrust = x(3*N+4:4*N+4);
final_time = x(4*N+5);

% Initial guess values for decision variables
x0(1:N+1) = 0;
x0(N+2:2*N+2) = 0;
x0(2*N+3:3*N+3) = m0;
x0(3*N+4:4*N+4) = m0*g0*2;
x0(4*N+5) = 0;

% % Initial guess values for decision variables
% x0(1:N+1) = 0;       % for altitude
% x0(N+2:2*N+2) = 0;    % for velocity
% x0(2*N+3:3*N+3) = m0;   % for mass
% x0(3*N+4:4*N+4) =  m0*g0*2;   % for thrust
% x0(4*N+5) = 0;             % for final time

% % Initial guess values for decision variables
% x0(1:N+1) = 3.1178e+04;       % for altitude
% x0(N+2:2*N+2) = 255.7725;    % for velocity
% x0(2*N+3:3*N+3) = 2.8203e+03;   % for mass
% x0(3*N+4:4*N+4) =  4.7044e+04;   % for thrust
% x0(4*N+5) = 187.6558;             % for final time

% linear inequality and equality constraints
A = [];
b = [];
Aeq = [];
beq = [];

% Lower and Upper bounds for the variables
lb(1:N+1) = 0;
lb(N+2:2*N+2) = 0;
lb(2*N+3) = m0;
lb(2*N+4:3*N+3) = m0-mp0;
lb(3*N+4:4*N+4) = -m0*g0*2;
lb(4*N+5) = 0;

ub(1:N+1) = inf;
ub(N+2:2*N+3) = inf;
ub(2*N+4:3*N+3) = m0 ;
ub(3*N+4:4*N+4) = m0*g0*2;
ub(4*N+5) = inf;

tic;
% optiSolver('NLP')
% Then select it via optiset:

opts = optiset('solver','IPOPT','maxiter',5000,'maxfeval',200000,'tolrfun',1e-7,'tolafun',1e-7, ...
            'display','iter');
% options =  optimoptions ('fmincon','Algorithm','sqp','Display','iter','OptimalityTolerance',...
% 1e-10 , 'ConstraintTolerance' ,1e-5, 'MaxIterations', 2000,'MaxFunctionEvaluations',...
% 200000);
   
    if strcmp(PS_method,'LGL')
       [x,fval,ef,output] = opti_fmincon(@(x) single_stage_objective_func(x,N),x0,A,b,Aeq,beq,lb,ub,@(x) single_stage_Nonlinear_func_LGL(x,N,D,problem),opts);
    elseif strcmp(PS_method,'LG') 
       [x,fval,ef,output] = fmincon(@(x) single_stage_objective_func(x,N),x0,A,b,Aeq,beq,lb,ub,@(x) single_stage_Nonlinear_func_LG(x,N,D,problem),options);
    elseif strcmp(PS_method,'LGR')
       [x,fval,ef,output] = fmincon(@(x) single_stage_objective_func(x,N),x0,A,b,Aeq,beq,lb,ub,@(x) single_stage_Nonlinear_func_LGR(x,N,D,problem),options);
    elseif strcmp(PS_method,'CGL')
       [x,fval,ef,output] = fmincon(@(x) single_stage_objective_func(x,N),x0,A,b,Aeq,beq,lb,ub,@(x) single_stage_Nonlinear_func_CGL(x,N,D,problem),options);
    end
    
% Stop the timer and display the elapsed time
elapsedTime = toc;
disp(['Elapsed time: ' num2str(elapsedTime) ' seconds']);

h = x(1:N+1);
v = x(N+2:2*N+2);
mass = x(2*N+3:3*N+3);
Thrust = x(3*N+4:4*N+4);
final_time = x(4*N+5);

% Lagrange interpolation
t = ((x(4*N+5)-t0)/2).*nodes+(x(4*N+5)+t0)/2;
z = 0:0.1:x(4*N+5);  % at time in seconds


collocation_points=t';
function_value=h;
altitude = lagrange_interpolation_n(collocation_points, function_value, z);

function_value=v;
velocity=lagrange_interpolation_n(collocation_points, function_value, z);

function_value=mass;
mass=lagrange_interpolation_n(collocation_points, function_value, z);

function_value = Thrust;
Thrust = lagrange_interpolation_n(collocation_points,function_value,z);

% figure
figure(1)
plot(z,altitude/1000,'g-','LineWidth',1.5)
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
plot(z,velocity/1000,'g-','LineWidth',1.5 )
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
plot(z,mass,'g-','LineWidth',1.5)
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
plot(z,Thrust/1000,'g-','LineWidth',1.5)
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
 

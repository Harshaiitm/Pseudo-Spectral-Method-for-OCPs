clc;clear all; close all;
%==============================================================================================%
%--- options ---%
% pseudospectral method
PS_method = 'LGL';   % either LGL or LG or LGR
N = 50;     % Order of the polynomial
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


x0(1:N) = 0;   % altitude
x0(N+1) = 19994.88; 
x0(N+2:2*N+1) = 129;     % velocity
x0(2*N+2) = 295.092;
x0(2*N+3:3*N+2) = 0;     % gamma
x0(3*N+3) = 0;          % final time
x0(3*N+4) = 324;
x0(3*N+5:4*N+5) = 19050.864;    % mass
x0(4*N+6:5*N+6) = -pi/4;      % alpha
     


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



% figure

figure(1)
plot(v,h,'k-')
ylim([0 20000])
xlabel('Airspeed [m/s]')
ylabel('Altitude [m]')
grid on

figure(2)
plot(t,gamma*180/pi,'k-')
xlim([0 tf])
xlabel('Time [s]')
ylabel('Flight Path Angle [deg]')
grid on

figure(3)
plot(t,h,'b-' )
xlim([0 tf])
ylim([0 20000])
xlabel('Time [s]')
ylabel('Altitude [m]')
grid on

figure(4)
plot(t,v,'b-')
xlabel('Time [s]')
ylabel('Velocity [m/s]')
grid on
 
figure(5)
plot(t,mass,'b-')
xlabel('Time [s]')
ylabel('Aircraft Mass [kg]')
grid on

figure(6)
plot(t,alpha,'b-' )
xlabel('Time [s]')
ylabel('Control Input (angle of attack) [deg]')
grid on



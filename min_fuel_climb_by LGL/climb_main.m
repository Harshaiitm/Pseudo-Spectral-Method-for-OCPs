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
tf = 500;
t = ((tf-t0)/2).*nodes+(tf+t0)/2;


x = zeros(5*N+5);
h = x(1:N+1);
v = x(N+2:2*N+2);
gamma = x(2*N+3:3*N+3);
mass = x(3*N+4:4*N+4);
alpha = x(4*N+5:5*N+5);
time = x(5*N+6);

r = h + Re;
[rho,sos]=atm_data(h);
Mach =  v./sos;


[Clalpha,CD0,eta] = aero_data(Mach);
Thrust = thrust_avialble(h,Mach);

CD = CD0 + eta.*(Clalpha.*alpha).^2;
CL = Clalpha.*alpha;
q = 0.5.*rho.*v.*v;
Drag = q.*S.*CD;
Lift = q.*S.*CL;


% Boundary Conditions 
h0 = 0; 
hf = 19994.88; 
v0 = 129.314; 
vf = 295.092;
gamma0 = 0; 
gammaf = 0; 
mass0 = 19050.864;


x0(1:N+1) = 0;         % altitude
x0(N+2:2*N+2) = 129;     % velocity
x0(2*N+3:3*N+3) = -40*pi/180;   % gamma
x0(3*N+4:4*N+4) = 22;         % mass
x0(4*N+5:5*N+5) = -pi/4;      % alpha
x0(5*N+6) = 324;


% linear inequality and equality constraints
A = [];
b = [];
Aeq = [];
beq = [];


% Lower and Upper bounds for the variables
lb(1:N+1) = 0;
lb(N+2:2*N+2) = 5;
lb(2*N+3:3*N+3) = -40*pi/180 ;
lb(3*N+4:4*N+4) = 22;
lb(4*N+5:5*N+5) = -pi/4;
lb(5*N+6) = 0;
ub(1:N+1) = 21031.2;
ub(N+2:2*N+2) = 1000;
ub(2*N+3:3*N+3) = 40*pi/180 ;
ub(3*N+4:4*N+4) = 20410;
ub(4*N+5:5*N+5) = pi/4;
ub(5*N+6) = 500;


%==============================================================================================%
% solver fmincon
% Start the timer
tic;

options =  optimoptions ('fmincon','Algorithm','sqp','Display','iter','OptimalityTolerance',...
1e-10 , 'ConstraintTolerance' ,1e-5, 'MaxIterations', 2000,'MaxFunctionEvaluations',...
20000);
   
[x,fval,ef,output] = fmincon(@(x) climb_objective_func(x,N),x0,A,b,Aeq,beq,lb,ub,@(x) climb_Nonlinear_func(x,N,D,mu,Thrust,r,t0,tf,Lift,Drag,g0,Isp),options);
   

% Stop the timer and display the elapsed time
elapsedTime = toc;
disp(['Elapsed time: ' num2str(elapsedTime) ' seconds']);

h = x(1:N+1);
v = x(N+2:2*N+2);
gamma = x(2*N+3:3*N+3);
mass = x(3*N+4:4*N+4);
alpha = x(4*N+5:5*N+5);
time = x(5*N+6);

%% figure

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
xlim([0 tf])?
ylim([0 20000])
xlabel('Time [s]')
ylabel('Altitude [m]')
grid on

% figure(4)
% plot(xx,speval(solution,'X',2,xx)/100,'b-' )
% plot(tv,xv(:,2)/100,'k-.')
% xlim([0 solution.tf])
% xlabel('Time [s]')
% ylabel('Velocity [100 m/s]')
% grid on
% 
% figure
% hold on
% plot(xx,speval(solution,'X',4,xx),'b-' )
% plot(tv,xv(:,4),'k-.')
% xlim([0 solution.tf])
% xlabel('Time [s]')
% ylabel('Aircraft Mass [kg]')
% grid on
% 
% figure
% hold on
% plot(xx,speval(solution,'U',1,xx),'b-' )
% plot(tv,uv(:,1),'k-.')
% xlim([0 solution.tf])
% xlabel('Time [s]')
% ylabel('Control Input (angle of attack) [deg]')
% grid on



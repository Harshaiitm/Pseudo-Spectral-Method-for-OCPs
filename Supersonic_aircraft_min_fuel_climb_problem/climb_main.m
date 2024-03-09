clc;clear all; close all;
%==============================================================================================%
%--- options ---%
% pseudospectral method
PS_method = 'LG';                          % either LGL or LG or LGR or CGL
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
        nodes = flip(nodes);                % Flipped LG method
        weights = flip(weights);            % weights are flipped
        nodes = [-1;nodes];                 % Introducing non-collocated point -1
        D=collocD(nodes);                   % segment differentiation matrix
        D(1,:) = [];                        % deletion of first row associated with non-collocated point
        nodes(1) = [];
    
    elseif  strcmp(PS_method,'CGL')
        N = M-1;                             % Order of the polynomial
        [nodes] = CGL_nodes(N);              % calculate scaled node locations and weights
        weights = CGL_weights(nodes);
        D=collocD(nodes);                    % segment differentiation matrix  
    end   
%================================================================================================================%
% Problem data    
Re = 6378145;
mu = 3.986e14;
S = 49.2386;
g0 = 9.80665;
Isp = 1600;
t0 = 0;
tf = 400;
t = ((tf-t0)/2).*nodes+(tf+t0)/2;

problem.Re = Re;
problem.mu = mu;
problem.S = S;
problem.g0 =g0;
problem.Isp = Isp;
problem.t0 = t0;
problem.weights = weights;


x = zeros(5*M+1);
h = x(1:M);
v = x(M+1:2*M);
gamma = x(2*M+1:3*M);
mass = x(3*M+1:4*M);
alpha = x(4*M+1:5*M);
final_time = x(5*M+1);

% Guess values
x0(1:M-1) = 0;                              % altitude
x0(M) = 19994.88;                           %  final altitude
x0(M+1:2*M-1) = 129.314;                    % initial velocity
x0(2*M) = 295.092;                          % final velocity
x0(2*M+1:3*M) = 0;                          %  gamma
x0(3*M+1:4*M) = 19050.864;                  % mass
x0(4*M+1:5*M) = 0;                          % alpha
x0(5*M+1) = 324;                            % final time


% linear inequality and equality constraints
A = [];
b = [];
Aeq = [];
beq = [];


% Lower and Upper bounds for the variables
lb(1:M) = 0;
lb(M+1:2*M) = 5;
lb(2*M+1:3*M) = -40*pi/180 ;
lb(3*M+1:4*M) = 22;
lb(4*M+1:5*M) = -pi/4;
lb(5*M+1) = 0;


ub(1:M) = 21031.2;
ub(M+1:2*M) = 1000;
ub(2*M+1:3*M) = 40*pi/180 ;
ub(3*M+1:4*M) = 20410;
ub(4*M+1:5*M) = pi/4;
ub(5*M+1) = 400;




%==============================================================================================%
% solver fmincon
% Start the timer
tic;

options =  optimoptions ('fmincon','Algorithm','sqp','Display','iter','OptimalityTolerance',...
1e-10 , 'ConstraintTolerance' ,1e-5, 'MaxIterations', 20000,'MaxFunctionEvaluations',...
200000);
   
    if strcmp(PS_method,'LGL')
       [x,fval,ef,output] = fmincon(@(x) climb_objective_func(x,M),x0,A,b,Aeq,beq,lb,ub,@(x) climb_Nonlinear_func_LGL(x,M,D,problem),options);
    elseif strcmp(PS_method,'LGR')
       [x,fval,ef,output] = fmincon(@(x) climb_objective_func(x,M),x0,A,b,Aeq,beq,lb,ub,@(x) climb_Nonlinear_func_LGR(x,x0,M,D,problem),options);
    elseif strcmp(PS_method,'LG') 
       [x,fval,ef,output] = fmincon(@(x) climb_objective_func(x,M),x0,A,b,Aeq,beq,lb,ub,@(x) climb_Nonlinear_func_LG(x,x0,M,D,problem),options);
    elseif strcmp(PS_method,'CGL')
       [x,fval,ef,output] = fmincon(@(x) climb_objective_func(x,M),x0,A,b,Aeq,beq,lb,ub,@(x) climb_Nonlinear_func_CGL(x,M,D,problem),options);
    end
   

% Stop the timer and display the elapsed time
elapsedTime = toc;
disp(['Elapsed time: ' num2str(elapsedTime) ' seconds']);

h = x(1:M);
v = x(M+1:2*M);
gamma = x(2*M+1:3*M);
mass = x(3*M+1:4*M);
alpha = x(4*M+1:5*M);
final_time = x(5*M+1);

if strcmp(PS_method,'LG')
    h = x(1:M);                               % altitude
    h(M+1)= x0(1) + weights'*h';
    v = x(M+1:2*M);                           % velocity  
    v(M+1)= x0(M+1) + weights'*v';
    gamma = x(2*M+1:3*M);                     % flight path angle
    gamma(M+1)= x0(2*M+1) + weights'*gamma';
    mass = x(3*M+1:4*M);                      % aircraft mass
    mass(M+1)= x0(3*M+1) + weights'*mass';
    alpha = x(4*M+1:5*M);                     % angle of attack
    final_time = x(5*M+1);
 end




% Lagrange interpolation
t = ((final_time-t0)/2).*nodes+(final_time+t0)/2;

z = t0:0.1:final_time;                  % at time in seconds

collocation_points = t';
function_value = h;
if strcmp(PS_method,'LGR')
    collocation_points = [t0 t'];
    function_value = [x0(1) h];
end
if strcmp(PS_method,'LG')
    collocation_points = [t0 t'];
    function_value = [x0(1) h];
end
altitude = lagrange_interpolation_n(collocation_points, function_value, z);


function_value = v;
if strcmp(PS_method,'LGR')
    function_value = [x0(M+1) v];
end
if strcmp(PS_method,'LG')
    function_value = [x0(M+1) v];
end
velocity = lagrange_interpolation_n(collocation_points, function_value, z);

function_value = gamma;
if strcmp(PS_method,'LGR')
    function_value= [x0(2*M+1) gamma];
end
if strcmp(PS_method,'LG')
    function_value= [x0(2*M+1) gamma];
end
gamma = lagrange_interpolation_n(collocation_points, function_value, z);

function_value = mass;
if strcmp(PS_method,'LGR')
    function_value = [x0(3*M+1) mass];
end
if strcmp(PS_method,'LG')
    function_value = [x0(3*M+1) mass];
end
mass = lagrange_interpolation_n(collocation_points,function_value,z);

function_value = alpha;
if strcmp(PS_method,'LGR')
    function_value= [x0(4*M+1) alpha];
end
if strcmp(PS_method,'LG')
    function_value= [x0(4*M+1) alpha];
end
alpha = lagrange_interpolation_n(collocation_points,function_value,z);


% figure
figure(1)
plot(velocity/100,altitude,'g-','LineWidth',1.5)
ylim([0 20000])
xlabel('Airspeed [100 m/s]')
ylabel('Altitude [m]')
hold on
load alt_VS_airspeed.csv
ai1 = alt_VS_airspeed(:,1);
ai2 = alt_VS_airspeed(:,2);
plot(ai1,ai2,'r--','LineWidth',2)
legend("PS Method","ICLOCS2");
% title("Altitude variation w.r.t airspeed",PS_method)
set(gca, 'FontSize', 40);
hold off
grid on

figure(2)
plot(z,altitude,'g-','LineWidth',2)
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
% title("Altitude variation w.r.t time",PS_method)
set(gca, 'FontSize', 40);
hold off 

figure(3)
plot(z,rad2deg(gamma),'g-','LineWidth',2)
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
% title("Flight Path angle variation w.r.t time",PS_method)
set(gca, 'FontSize', 40);
hold off



figure(4)
plot(z,velocity/100,'g-','LineWidth',2)
xlabel('Time [s]')
ylabel('Velocity [100 m/s]')
grid on
hold on
load velocity_vs_time.csv
al1 = velocity_vs_time(:,1);
al2 = velocity_vs_time(:,2);
plot(al1,al2,'r--','LineWidth',1.5);
grid on
legend("PS Method","ICLOCS2");
% title("Aircraft velocity variation w.r.t time",PS_method)
set(gca, 'FontSize', 40);
hold off 
 
figure(5)
plot(z,mass,'g-','LineWidth',2)
xlabel('Time [s]')
ylabel('Aircraft Mass [kg]')
grid on
hold on
load mass_vs_time.csv
al1 = mass_vs_time(:,1);
al2 = mass_vs_time(:,2);
plot(al1,al2,'r--','LineWidth',1.5);
legend("PS Method","ICLOCS2");
% title("Aircraft mass variation w.r.t time",PS_method)
set(gca, 'FontSize', 40);
grid on
hold off 

figure(6)
plot(z,rad2deg(alpha),'g-','LineWidth',2 )
xlabel('Time [s]')
ylabel('angle of attack [deg]')
grid on
hold on
load alpha_vs_time.csv
al1 = alpha_vs_time(:,1);
al2 = alpha_vs_time(:,2);
plot(al1,rad2deg(al2),'r--','LineWidth',1.5);
legend("PS Method","ICLOCS2");
% title("Control input(alpha) variation w.r.t time",PS_method)
set(gca, 'FontSize', 40);
hold off



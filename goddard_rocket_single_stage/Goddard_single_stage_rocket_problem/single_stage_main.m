% single stage Goddard Rocket Problem
% single_stage_main.m
clc;clear all; close all;
%==============================================================================================%
%--- options ---%
% pseudospectral method
PS_method = 'LGL';                           % either LGL or LG or LGR or CGL
M = 45;                                     % Number of collocation points
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


% Decision veriables
x = zeros(4*M+1);
h = x(1:M);
v = x(M+1:2*M);
mass = x(2*M+1:3*M);
Thrust = x(3*M+1:4*M);
final_time = x(4*M+1);

% Initial guess values for decision variables
x0(1:M) = 0;
x0(M+1:2*M) = 0;
x0(2*M+1:3*M) = m0;
x0(3*M+1:4*M) = m0*g0*2;
x0(4*M+1) = 20;

% Initial and Final Conditions
hi = 0;
vi = 0;
mass_i = m0;
mass_f = m0 - mp0;
Thrust_i = m0*g0*2;

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
problem.hi = hi;
problem.vi = vi;
problem.mass_i = mass_i;
problem.mass_f = mass_f;
problem.Thrust_i = Thrust_i;
problem.x0 = x0;



% linear inequality and equality constraints
A = [];
b = [];
Aeq = [];
beq = [];

% Lower and Upper bounds for the variables
lb(1:M) = 0;
lb(M+1:2*M) = 0;
lb(2*M+1) = m0;
lb(2*M+2:3*M) = m0-mp0;
lb(3*M+1:4*M) = 0;
lb(4*M+1) = 0;

ub(1:M) = inf;
ub(M+1:2*M) = inf;
ub(2*M+1:3*M) = m0 ;
ub(3*M+1:4*M) = m0*g0*2;
ub(4*M+1) = inf;

tic;
% optiSolver('NLP')
% Then select it via optiset:

% opts = optiset('solver','IPOPT','maxiter',5000,'maxfeval',200000,'tolrfun',1e-7,'tolafun',1e-7, ...
%             'display','iter');
options =  optimoptions ('fmincon','Algorithm','sqp','Display','iter','OptimalityTolerance',...
1e-10 , 'ConstraintTolerance' ,1e-5, 'MaxIterations', 2000,'MaxFunctionEvaluations',...
200000);
   
    if strcmp(PS_method,'LGL')
       [x,fval,ef,output] = fmincon(@(x) single_stage_objective_func(x,M),x0,A,b,Aeq,beq,lb,ub,@(x) single_stage_Nonlinear_func_LGL(x,M,D,problem),options);
    elseif strcmp(PS_method,'LGR')
       [x,fval,ef,output] = fmincon(@(x) single_stage_objective_func(x,M),x0,A,b,Aeq,beq,lb,ub,@(x) single_stage_Nonlinear_func_LGR(x,M,D,problem),options);
    elseif strcmp(PS_method,'LG') 
       [x,fval,ef,output] = fmincon(@(x) single_stage_objective_func(x,M),x0,A,b,Aeq,beq,lb,ub,@(x) single_stage_Nonlinear_func_LG(x,M,D,problem),options);
    elseif strcmp(PS_method,'CGL')
       [x,fval,ef,output] = fmincon(@(x) single_stage_objective_func(x,M),x0,A,b,Aeq,beq,lb,ub,@(x) single_stage_Nonlinear_func_CGL(x,M,D,problem),options);
    end
    
% Stop the timer and display the elapsed time
elapsedTime = toc;
disp(['Elapsed time: ' num2str(elapsedTime) ' seconds']);

h = x(1:M);
v = x(M+1:2*M);
mass = x(2*M+1:3*M);
Thrust = x(3*M+1:4*M);
final_time = x(4*M+1);

if strcmp(PS_method,'LG')
    h = x(1:M);
    h(M+1) = hi + weights'*h';
    v = x(M+1:2*M);
    v(M+1) = vi + weights'*v';
    mass = x(2*M+1:3*M);
    mass(M+1) = mass_i + weights'*mass';
    Thrust = x(3*M+1:4*M);
    final_time = x(4*M+1); 
end

% Lagrange interpolation
t = ((final_time-t0)/2).*nodes+(final_time+t0)/2;
z = t0:0.1:final_time;  % at time in seconds

collocation_points = t';
function_value = h;

if strcmp(PS_method,'LGR')
    collocation_points=[t0 t'];
    function_value= [hi h];
end
if strcmp(PS_method,'LG')
    collocation_points=[t0 t'];
    function_value= [hi h];
end
altitude = lagrange_interpolation_n(collocation_points, function_value, z);

function_value=v;
if strcmp(PS_method,'LGR')
    function_value= [vi v];
end
if strcmp(PS_method,'LG')
    collocation_points=[t0 t'];
    function_value= [vi v];
end
velocity=lagrange_interpolation_n(collocation_points, function_value, z);


function_value=mass;
if strcmp(PS_method,'LGR')
    collocation_points=[t0 t'];
    function_value= [mass_i mass];
end
if strcmp(PS_method,'LG')
    collocation_points=[t0 t'];
    function_value= [mass_i mass];
end
mass=lagrange_interpolation_n(collocation_points, function_value, z);

if strcmp(PS_method,'LGR')
   collocation_points = t'; 
end
if strcmp(PS_method,'LG')
    collocation_points= t';
end
function_value = Thrust;
Thrust = lagrange_interpolation_n(collocation_points,function_value,z);

% figure
figure(1)
plot(z,altitude/1000,'k-','LineWidth',1.5)
xlabel('Time [s]')
ylabel('Altitude [km]')
hold on
load alt_VS_time.csv
ai1 = alt_VS_time(:,1);
ai2 = alt_VS_time(:,2);
plot(ai1,ai2,'r--','LineWidth',2)
legend("PS Method","NPSOL");
title("Altitude variation w.r.t time",PS_method)
set(gca, 'FontSize', 20);
hold off
grid on

figure(2)
plot(z,velocity/1000,'k-','LineWidth',1.5)
xlabel('Time [s]')
ylabel('velocity [km/s]')
hold on
load velocity_vs_time.csv
al1 = velocity_vs_time(:,1);
al2 = velocity_vs_time(:,2);
plot(al1,al2,'r--','LineWidth',1.5);
grid on
legend("PS Method","NPSOL");
title("Velocity variation w.r.t time",PS_method)
set(gca, 'FontSize', 20);
hold off 

figure(3)
plot(z,mass,'k-','LineWidth',2)
xlabel('Time [s]')
ylabel('Mass [kg]')
grid on
hold on
load mass_vs_time.csv
al1 = mass_vs_time(:,1);
al2 = mass_vs_time(:,2);
plot(al1,al2,'r--','LineWidth',2);
grid on
legend("PS Method","NPSOL");
title("Vehicle mass variation w.r.t time",PS_method)
set(gca, 'FontSize', 20);
hold off



figure(4)
plot(z,Thrust/1000,'k-','LineWidth',1.5)
ylim([-2,100])
xlabel('Time [s]')
ylabel('Thrust [kN]')
grid on
hold on
load Thrust_vs_time.csv
al1 = Thrust_vs_time(:,1);
al2 = Thrust_vs_time(:,2);
plot(al1,al2,'r--','LineWidth',2);
grid on
legend("PS Method","NPSOL");
title("Thrust variation w.r.t time",PS_method)
set(gca, 'FontSize', 20);
hold off 
 
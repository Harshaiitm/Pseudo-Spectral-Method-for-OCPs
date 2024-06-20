%--------------------------------------------------------------------------
% DImain.m
% - Main script for the Double integrator Tracking problem
% - Define the Pseudo spectral method and order of the polynomial then solve and plot solution
%--------------------------------------------------------------------------
% DI_main
%--------------------------------------------------------------------------
clc; clear all; close all;
%==============================================================================================%
%--- options ---%
% pseudospectral method
PS_method = 'CGL';                           % either LGL or LG or LGR or CGL
M = 9;                                     % Number of collocation points 
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
        nodes = [-1;nodes];               % Introducing non-collocated point -1 and 1
        D=collocD(nodes);                   % segment differentiation matrix
        D(1,:) = [];                        % deletion of first row associated with non-collocated point
        nodes(1) = [];

    elseif  strcmp(PS_method,'CGL')
        N = M-1;                             % Order of the polynomial
        [nodes] = CGL_nodes(N);              % calculate scaled node locations and weights
        weights = CGL_weights(nodes);
        D=collocD(nodes);                    % segment differentiation matrix  
    end   
    
%==============================================================================================%
x = zeros(1,3*M);                            % state and control vector assigned       
x1 = x(1:M);                                 % position                        
x2 = x(M+1:2*M);                             % velocity  
x3 = x(2*M+1:3*M);                           % accleration


% initialization of state and control variables
t0 = 0;                                     % initial time   
tf = 10;                                    % final time
t = ((tf-t0)/2).*nodes+(tf+t0)/2;

x0(1:M) = 0;                                % position
x0(M+1:2*M) = 5;                            % velocity
x0(2*M+1:3*M) = 0;                          % (Force for unit mass)

% starting and ending points
xi = 0;                                     % position                            
vi = 5;                                     % velocity

problem.xi = xi;
problem.vi = vi;
problem.x0 = x0;
problem.t0 = t0;
problem.tf = tf;
problem.weights = weights;


% linear inequality and equality constraints
A = [];
b = [];
Aeq = [];
beq = [];

% Lower and Upper bounds for the variables
lb(1:M) = -6;
lb(M+1:2*M) = -10;
lb(2*M+1:3*M) = -10;
ub(1:M) = 6;
ub(M+1:2*M) = 10;
ub(2*M+1:3*M) = 10;



%==============================================================================================%
% solver fmincon
% Start the timer
tic;

options =  optimoptions ('fmincon','Algorithm','sqp','Display','iter','OptimalityTolerance',...
1e-10 , 'ConstraintTolerance' ,1e-5, 'MaxIterations', 20000,'MaxFunctionEvaluations',...
500000);
    if strcmp(PS_method,'LGL')
       [x,fval,ef,output] = fmincon(@(x) DI_Objective_func(x,M,weights,t,problem),x0,A,b,Aeq,beq,lb,ub,@(x) DI_Nonlinearcon_LGL(x,M,D,problem),options);
    elseif strcmp(PS_method,'LGR')
       [x,fval,ef,output] = fmincon(@(x) DI_Objective_func(x,M,weights,t,problem),x0,A,b,Aeq,beq,lb,ub,@(x) DI_Nonlinearcon_LGR(x,M,D,problem),options);
    elseif strcmp(PS_method,'LG') 
       [x,fval,ef,output] = fmincon(@(x) DI_Objective_func(x,M,weights,t,problem),x0,A,b,Aeq,beq,lb,ub,@(x) DI_Nonlinearcon_LG(x,M,D,problem),options);
    elseif strcmp(PS_method,'CGL')
       [x,fval,ef,output] = fmincon(@(x) DI_Objective_func(x,M,weights,t,problem),x0,A,b,Aeq,beq,lb,ub,@(x) DI_Nonlinearcon_CGL(x,M,D,problem),options);
    end 

% Stop the timer and display the elapsed time
elapsedTime = toc;
disp(['Elapsed time: ' num2str(elapsedTime) ' seconds']);



%========================================================================================================
% Plotting                                     % state and control vector assigned       
x1 = x(1:M);                                   % position                        
x2 = x(M+1:2*M);                               % velocity  
x3 = x(2*M+1:3*M);                             % accleration

if strcmp(PS_method,'LG')
    x1 = x(1:M);
    x2 = x(M+1:2*M);
    x3 = x(2*M+1:3*M); 
    xf = 5*sin(tf); 
    vf = 5*cos(tf);
end
vf = 5*cos(tf);

% Lagrange interpolation
collocation_points=t';

tic;
z = t0:0.01:tf;
function_value= x1;

if strcmp(PS_method,'LGR')
    collocation_points=[t0 t'];
    function_value= [xi x1];
end
if strcmp(PS_method,'LG')
    collocation_points=[t0 t' tf];
    function_value= [xi x1 xf];
end
x1R = 5*sin(z);
x2R = 5*cos(z);
position= lagrange_interpolation_n(collocation_points, function_value, z);

function_value=x2;
if strcmp(PS_method,'LGR')
    function_value= [vi x2 vf];
end
if strcmp(PS_method,'LG')
    function_value= [vi x2 vf];
end
velocity=lagrange_interpolation_n(collocation_points, function_value, z);

if strcmp(PS_method,'LGR')
   collocation_points = t'; 
end
if strcmp(PS_method,'LG')
   collocation_points = t';
end
function_value=x3;
acceleration=lagrange_interpolation_n(collocation_points, function_value, z);

elapsedTime = toc;
disp(['Elapsed time: ' num2str(elapsedTime) ' seconds']);


figure(1)
plot(z, position,'LineWidth', 2.5);
hold on
plot(z, x1R,'LineWidth', 2.5);
xlabel('Time (s)');
ylabel('Position (m)');
% title('Double Integrator Tracking Problem',PS_method);
% legend({'actual position','reference position'},Location="northeast");
set(gca, 'FontSize', 40);
grid on
hold off

figure(2)
plot(z, velocity,'LineWidth', 2.5);
hold on
plot(z, x2R,'LineWidth', 2.5);
xlabel('Time (s)');
ylabel('Velocity (m/s)');
% title('Double Integrator Tracking Problem',PS_method);
% legend({'Actual','Reference'},Location="northeast");
set(gca, 'FontSize', 40);
grid on
hold off

figure(3)
plot(z,acceleration,'LineWidth', 2.5);
xlabel('Time (s)');
ylabel('control variable (N)');
% legend({'control variable'},Location="northeast");
title('Double integrator tracking problem',PS_method);
grid on
set(gca, 'FontSize', 40);









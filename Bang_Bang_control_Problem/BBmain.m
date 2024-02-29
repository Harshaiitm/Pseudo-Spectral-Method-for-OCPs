%--------------------------------------------------------------------------
% BBmain.m
% - Main script for the Bang_Bang control problem
% - Define the PS method and Order of the polynomial, then solve and plot solution
%--------------------------------------------------------------------------
% BBmain
%--------------------------------------------------------------------------
clc; clear all; close all;
%==============================================================================================%
%--- options ---%
% pseudospectral method
PS_method = 'LGR';   % either LGL or LG or LGR or CGL
N = 100;     % Order of the polynomial
addpath('.../PS_methods')  % add the PS_method file directory

     if  strcmp(PS_method,'LGL')
        [nodes,weights] = LGL_nodes(N); % calculate scaled node locations and weights
        D=collocD(nodes); % segment differentiation matrix
    elseif strcmp(PS_method,'LG')
        nodes(1) =-1;
        [nodes(2:N+1,1),weights]=LG_nodes(N,-1,1); % calculate scaled node locations and weights
        nodes(2:N+1,1) = flip(nodes(2:N+1,1)); 
        D=collocD(nodes); % segment differentiation matrix
        D(N+1,:) = [];
    elseif strcmp(PS_method,'LGR')
        [nodes,weights] = LGR_nodes(N); % calculate scaled node locations and weights
        nodes = flip(-nodes);
        nodes(1) = -1;
        D = collocD(nodes); % segment differentiation matrix
        weights = flip(weights);
        D(N+1,:) = [];
    elseif  strcmp(PS_method,'CGL')
        [nodes] = CGL_nodes(N);     % calculate scaled node locations and weights
        weights = CGL_weights(nodes);
        D=collocD(nodes);           % segment differentiation matrix  
    end       
   

%==============================================================================================%

x = zeros(1,3*N+4); 
x1 = x(1:N+1);              % position vector
x2 = x(N+2:2*N+2);          % velocity vector
x3 = x(2*N+3:3*N+3);        % Acceleration vector
x4 = x(3*N+4);              % final time
t0 = 0;
tf = 35;
t = ((tf-t0)/2).*nodes+(tf+t0)/2;

% initialization of state and control variables
x0(1) = 0;
x0(2:N+1) = 0;
x0(N+2) = 0;
x0(N+3:2*N+1) = 0;
x0(2*N+2) = 0;
x0(2*N+3:3*N+3) =1;
x0(3*N+4) = 5;


A = [];
b = [];
Aeq = [];
beq = [];
lb = x;
lb(1:N+1) = -10;
lb(N+2:2*N+2) = -200;
lb(2*N+3:3*N+3) = -2;
lb(3*N+4) = 0;
ub = x;
ub(1:N+1) = 300;
ub(N+2:2*N+2) = 200;
ub(2*N+3:3*N+3) = 1;
ub(3*N+4) = 35;


% Start the timer
tic;

options =  optimoptions ('fmincon','Algorithm','sqp','Display','iter','OptimalityTolerance',...
1e-10 , 'ConstraintTolerance' ,1e-5, 'MaxIterations', 2000,'MaxFunctionEvaluations',...
200000);
    if strcmp(PS_method,'LGL')
       [x,fval,ef,output] = fmincon(@(x) Objective_func(x,N),x0,A,b,Aeq,beq,lb,ub,@(x) Nonlinearcon_LGL(x,N,D,t0),options);
    elseif strcmp(PS_method,'LG') 
       [x,fval,ef,output] = fmincon(@(x) Objective_func(x,N),x0,A,b,Aeq,beq,lb,ub,@(x) Nonlinearcon_LG(x,N,D,t0,tf,weights),options);
    elseif strcmp(PS_method,'LGR')
       [x,fval,ef,output] = fmincon(@(x) Objective_func(x,N),x0,A,b,Aeq,beq,lb,ub,@(x) Nonlinearcon_LGR(x,N,D,t0,tf),options);
    elseif strcmp(PS_method,'CGL')
       [x,fval,ef,output] = fmincon(@(tf) Objective_func(x,N),x0,A,b,Aeq,beq,lb,ub,@(x) Nonlinearcon_CGL(x,N,D,t0,tf),options);
    end 

% Stop the timer and display the elapsed time
elapsedTime = toc;
disp(['Elapsed time: ' num2str(elapsedTime) ' seconds']);



%========================================================================================================%
% Plots

x1 = x(1:N+1);
x2 = x(N+2:2*N+2);
x3 = x(2*N+3:3*N+3);
x4 = x(3*N+4);
t = ((x4-t0)/2).*nodes+(x4+t0)/2;


% Lagrange interpolation
z = 0:0.1:x4;

collocation_points = t;
function_value=x1;
position= lagrange_interpolation_n(collocation_points, function_value, z);


function_value=x2;
velocity=lagrange_interpolation_n(collocation_points, function_value, z);


function_value=x3;
acceleration=lagrange_interpolation_n(collocation_points, function_value, z);




figure(1)
plot(z, position,'LineWidth',2);
hold on
plot(z, velocity,'LineWidth',2);
xlabel('Time (s)');
ylabel('State Variables');
xlim([0,30])
title('Final time minimization problem',PS_method);
legend({'Positon(x1)','Velocity(x2)'},Location="northeast");
% set(gca, 'FontSize', 40);
grid on
hold off


figure(2)
plot(z,acceleration,'LineWidth',1);
xlabel('Time (s)');
ylabel('countrol variables (N)');
legend({'control variable'},Location="northeast");
title('Bang Bang control problem',PS_method);
xlim([0,30])
grid on
% set(gca, 'FontSize', 40);




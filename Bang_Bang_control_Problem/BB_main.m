%--------------------------------------------------------------------------
% BB_main.m
% - Main script for the Bang_Bang control problem
% - Define the PS method and Order of the polynomial, then solve and plot solution
%--------------------------------------------------------------------------
% BB_main
%--------------------------------------------------------------------------
clc; clear all; close all;
%==============================================================================================%
%--- options ---%
% pseudospectral method
PS_method = 'CGL';                          % either LGL or LG or LGR or CGL
M = 100;                                     % number of collocation points

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
        nodes = [-1;nodes;1];               % Introducing non-collocated point -1
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
%==============================================================================================%

x = zeros(1,3*M+1); 
x1 = x(1:M);                        % position vector
x2 = x(M+1:2*M);                    % velocity vector
x3 = x(2*M+1:3*M);                  % Acceleration vector
x4 = x(3*M+1);                      % final time
t0 = 0;
% tf = 35;
% t = ((tf-t0)/2).*nodes+(tf+t0)/2;

% initialization of state and control variables
x0(1:M) = 0;
x0(M+1:2*M) = 0;
x0(2*M+1:3*M) = 1;
x0(3*M+1) = 35;


xi = 0;
xf = 300;
vi = 0;
vf = 0;

problem.xi = xi;
problem.xf = xf;
problem.vi = vi;
problem.vf = vf;
problem.x0 = x0;
problem.t0 = t0;

A = [];
b = [];
Aeq = [];
beq = [];

lb(1:M) = -10;
lb(M+1:2*M) = -200;
lb(2*M+1:3*M) = -2;
lb(3*M+1) = 0;

ub(1:M) = 300;
ub(M+1:2*M) = 200;
ub(2*M+1:3*M) = 1;
ub(3*M+1) = 35;


% Start the timer
tic;

options =  optimoptions ('fmincon','Algorithm','sqp','Display','iter','OptimalityTolerance',...
1e-10 , 'ConstraintTolerance' ,1e-5, 'MaxIterations', 2000,'MaxFunctionEvaluations',...
200000);
    if strcmp(PS_method,'LGL')
       [x,fval,ef,output] = fmincon(@(x) BB_Objective_func(x,M),x0,A,b,Aeq,beq,lb,ub,@(x) BB_Nonlinearcon_LGL(x,M,D,problem),options);
    elseif strcmp(PS_method,'LGR')
       [x,fval,ef,output] = fmincon(@(x) BB_Objective_func(x,M),x0,A,b,Aeq,beq,lb,ub,@(x) BB_Nonlinearcon_LGR(x,M,D,problem),options);
    elseif strcmp(PS_method,'LG') 
       [x,fval,ef,output] = fmincon(@(x) BB_Objective_func(x,M),x0,A,b,Aeq,beq,lb,ub,@(x) BB_Nonlinearcon_LG(x,M,D,problem),options);
    elseif strcmp(PS_method,'CGL')
       [x,fval,ef,output] = fmincon(@(x) BB_Objective_func(x,M),x0,A,b,Aeq,beq,lb,ub,@(x) BB_Nonlinearcon_CGL(x,M,D,problem),options);
    end 

% Stop the timer and display the elapsed time
elapsedTime = toc;
disp(['Elapsed time: ' num2str(elapsedTime) ' seconds']);



%========================================================================================================%
% Plots

x1 = x(1:M);
x2 = x(M+1:2*M);
x3 = x(2*M+1:3*M);
x4 = x(3*M+1);
t = ((x4-t0)/2).*nodes+(x4+t0)/2;

if strcmp(PS_method,'LG')
    x1 = x(1:M);                               % position                        
    x2 = x(M+1:2*M);                           % velocity  
    x3 = x(2*M+1:3*M);                         % accleration
    x4 = x(3*M+1);
end




% Lagrange interpolation
z = t0:0.1:x4;

collocation_points = t;
function_value = x1;


if strcmp(PS_method,'LGR')
    collocation_points = [t0 t'];
    function_value = [xi x1];
end
if strcmp(PS_method,'LG')
    collocation_points = [t0 t' x4];
    function_value = [xi x1 xf];
end
position = lagrange_interpolation_n(collocation_points, function_value, z);


function_value = x2;
if strcmp(PS_method,'LGR')
    collocation_points = [t0 t'];
    function_value = [vi x2];
end
if strcmp(PS_method,'LG')
    collocation_points = [t0 t' x4];
    function_value = [vi x2 vf];
end
velocity = lagrange_interpolation_n(collocation_points, function_value, z);


if strcmp(PS_method,'LGR')
   collocation_points = t'; 
end
if strcmp(PS_method,'LG')
   collocation_points = t';
end
function_value = x3;
acceleration = lagrange_interpolation_n(collocation_points, function_value, z);




figure(1)
plot(z, position,'LineWidth',1.5);
hold on
plot(z, velocity,'LineWidth',1.5);
xlabel('Time (s)');
ylabel('State Variables');
xlim([0,30])
title('Double Integrator Min-Time Repositioning',PS_method);
legend({'Positon(m)','Velocity(m/s)'},Location="northeast");
set(gca, 'FontSize', 40);
grid on
hold off


figure(2)
plot(z,acceleration,'LineWidth',1.5);
xlabel('Time (s)');
ylabel('countrol variable (N)');
legend({'control variable'},Location="northeast");
title('Double Integrator Min-Time Repositioning',PS_method);
xlim([0,30])
grid on
set(gca, 'FontSize', 40);




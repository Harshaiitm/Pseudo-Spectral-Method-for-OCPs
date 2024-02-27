%--------------------------------------------------------------------------
% DImain.m
% - Main script for the Double integrator Tracking problem
% - Define the Pseudo spectral method and order of the polynomial then solve and plot solution
%--------------------------------------------------------------------------
% DImain
%--------------------------------------------------------------------------
clc; clear all; close all;
%==============================================================================================%
%--- options ---%
% pseudospectral method
PS_method = 'LGR';   % either LGL or LG or LGR
N = 8;     % Order of the polynomial
addpath('C:\Users\Harshad\OneDrive\Desktop\Optimal Trajectory Transfer\DIT\PS_methods') % add the PS_method file directory

    if  strcmp(PS_method,'LGL')
        [nodes,weights,D_LGL] = LGL_computations(N); % calculate scaled node locations and weights
        % D=collocD(nodes); % segment differentiation matrix
        D = D_LGL;
    elseif strcmp(PS_method,'LG')
        nodes(1) =-1;
        [nodes(2:N+1,1),weights]=LG_nodes(N,-1,1); % calculate scaled node locations and weights
        nodes(2:N+1,1) = flip(nodes(2:N+1,1));
        nodes(N+1) = 1;
        D=collocD(nodes); % segment differentiation matrix
        D(N+1,:) = [];
    elseif strcmp(PS_method,'LGR')
        [nodes,weights] = LGR_nodes(N-1); % calculate scaled node locations and weights
        nodes = flip(-nodes);
        nodes = [-1;nodes];
        weights = flip(weights);
        D = collocD(nodes);          % segment differentiation matrix
        D(1,:) = [];
        % nodes(1) = []; 
    elseif  strcmp(PS_method,'CGL')
        [nodes] = CGL_nodes(N);     % calculate scaled node locations and weights
        weights = CGL_weights(nodes);
        D=collocD(nodes);           % segment differentiation matrix  
    end   
    
%==============================================================================================%
x = zeros(1,3*N+3);              % state and control vector assigned       
x1 = x(1:N+1);                   % position                        
x2 = x(N+2:2*N+2);               % velocity  
x3 = x(2*N+3:3*N+3);             % accleration
if strcmp(PS_method,'LG')
    x1(N+1)=x1(1)+weights(1:N)'*(x1(1:N))';
    x2(N+1)=x2(1)+weights(1:N)'*(x2(1:N))';
    x3(N+1)=x3(1)+weights(1:N)'*(x3(1:N))';
end

% initialization of state and control variables
t0 = 0;                          % initial time   
tf = 10;                         % final time
t = ((tf-t0)/2).*nodes+(tf+t0)/2;
x0(1:N+1) = 0;         % position
x0(N+2:2*N+2) = 5;     % velocity
x0(2*N+3:3*N+3) = 0;          % (Force for unit mass)


% linear inequality and equality constraints
A = [];
b = [];
Aeq = [];
beq = [];

% Lower and Upper bounds for the variables
lb(1:N+1) = -6;
lb(N+2:2*N+2) = -10;
lb(2*N+3:3*N+3) = -10;
ub(1:N+1) = 6;
ub(N+2:2*N+2) = 10;
ub(2*N+3:3*N+3) = 10;


%==============================================================================================%
% solver fmincon
% Start the timer
tic;

options =  optimoptions ('fmincon','Algorithm','sqp','Display','iter','OptimalityTolerance',...
1e-10 , 'ConstraintTolerance' ,1e-5, 'MaxIterations', 20000,'MaxFunctionEvaluations',...
500000);
    if strcmp(PS_method,'LGL')
       [x,fval,ef,output] = fmincon(@(x) Objective_LGL(x,N,weights,t0,tf,t),x0,A,b,Aeq,beq,lb,ub,@(x) Nonlinearcon_LGL(x,N,D,t0,tf),options);
    elseif strcmp(PS_method,'LG') 
       [x,fval,ef,output] = fmincon(@(x) Objective_LG(x,N,weights,t0,tf,t),x0,A,b,Aeq,beq,lb,ub,@(x) Nonlinearcon_LG(x,N,D,t0,tf,weights),options);
    elseif strcmp(PS_method,'LGR')
       [x,fval,ef,output] = fmincon(@(x) Objective_LGR(x,x0,N,weights,t0,tf,t),x0,A,b,Aeq,beq,lb,ub,@(x) Nonlinearcon_LGR(x,x0,N,D,t0,tf),options);
    elseif strcmp(PS_method,'CGL')
       [x,fval,ef,output] = fmincon(@(x) Objective_CGL(x,N,weights,t0,tf,t),x0,A,b,Aeq,beq,lb,ub,@(x) Nonlinearcon_CGL(x,N,D,t0,tf),options);
    end 

% Stop the timer and display the elapsed time
elapsedTime = toc;
disp(['Elapsed time: ' num2str(elapsedTime) ' seconds']);



%========================================================================================================
% Plotting               % state and control vector assigned       
x1 = x(1:N+1);                   % position                        
x2 = x(N+2:2*N+2);               % velocity  
x3 = x(2*N+3:3*N+3);             % accleration
if strcmp(PS_method,'LGR')
    x1 = [x0(1) x(2:N+1)];                       % position
    x2 = [x0(N+2) x(N+3:2*N+2)];                 % velocity
    x3 = [x(2*N+4:3*N+3)];             % accleration
end    
if strcmp(PS_method,'LG')
    x1 = x(1:N);                   % position                        
    x2 = x(N+2:2*N+1);               % velocity  
    x3 = x(2*N+2:3*N+2);             % accleration
    x1(N+1)=x1(1)+weights(1:N)'*(x1(1:N))';
    x2(N+1)=x2(1)+weights(1:N)'*(x2(1:N))';
    x3(N+1)=x3(1)+weights(1:N)'*(x3(1:N))';
end

% Lagrange interpolation
% z = linspace(1,1000,length(nodes));  % at time in seconds
collocation_points=t';
% z_value =t;


tic;
z = 0:0.01:t(end);
function_value=x1;
x1R = 5*sin(z);
x2R = 5*cos(z);
position= lagrange_interpolation_n(collocation_points, function_value, z);
% disp(['Position Equation:', char(Position_equation)]);
% position = subs(Position_equation,z,z_value);
% disp(['Position =',char(position),'m']);

function_value=x2;
velocity=lagrange_interpolation_n(collocation_points, function_value, z);
% % disp(['Velocity Equation:', char(velocity_equation)]);
% velocity = subs(velocity_equation,z,z_value);
% % disp(['Velocity =',char(velocity),'m/s']);

if strcmp(PS_method,'LGR')
   collocation_points = (t(2:N+1))'; 
end
function_value=x3;
acceleration=lagrange_interpolation_n(collocation_points, function_value, z);
% % disp(['Acceleration =',char(acceleration_equation),'m/s^2']);
% acceleration = subs(acceleration_equation,z,z_value);
% % disp(['acceleration =',char(acceleration),'m/s^2']);
% 
elapsedTime = toc;
disp(['Elapsed time: ' num2str(elapsedTime) ' seconds']);


figure(1)
plot(z, position,'LineWidth', 1.5);
hold on
plot(z, x1R,'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Position (m)');
title('Double Integrator Tracking Problem',PS_method);
legend({'actual position','reference position'},Location="northeast");
set(gca, 'FontSize', 40);
grid on
hold off

figure(2)
plot(z, velocity,'LineWidth', 1.5);
hold on
plot(z, x2R,'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Velocity (m/s)');
title('Double Integrator Tracking Problem',PS_method);
legend({'actual velocity','reference velocity'},Location="northeast");
set(gca, 'FontSize', 40);
grid on
hold off



figure(3)
plot(z,acceleration,'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('countrol variable (N)');
legend({'control variable'},Location="northeast");
title('Double integrator tracking problem',PS_method);
grid on
set(gca, 'FontSize', 40);










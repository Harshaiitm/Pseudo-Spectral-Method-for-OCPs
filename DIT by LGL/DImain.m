clc; clear all; close all;
%==============================================================================================%
N = 30; % Order of the polynomial
[nodes,weights] = LGL_nodes(N);  % Legendre_gauss_Lobatto PS Method 
D = collocD(nodes);              % Differential matrix for LGL                               
%==============================================================================================%
x = zeros(1,3*N+3);              % state and control vector assigned       
x1 = x(1:N+1);                   % position                        
x2 = x(N+2:2*N+2);               % velocity  
x3 = x(2*N+3:3*N+3);             % accleration


% initialization of state and control variables
t0 = 0;                          % initial time   
tf = 10;                         % final time
t = ((tf-t0)/2).*nodes+(tf+t0)/2;
x1 = 5*sin(t);
x2 = 5*cos(t);
% i = linspace(1,10,N+1)
% x1 = 5*sin(i);
% x2 = 5*cos(i);
x0(1:N+1) = double(x1);
x0(N+2:2*N+2) = double(x2);
x0(2*N+3:3*N+3) = 0.1;
 

% call for Objective and Nonlinearconstraints function
Objective = Objectivee_func(x,N,weights,t);
[c, ceq, dc, dceq] = Nonlinearcon(x,N,D,t0,tf);

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
1e-10 , 'ConstraintTolerance' ,1e-5, 'MaxIterations', 2000,'MaxFunctionEvaluations',...
20000);
[x,fval,ef,output] = fmincon(Objective,x0,A,b,Aeq,beq,lb,ub,@(x) Nonlinearcon(x,N,D,t0,tf),options);

% Stop the timer and display the elapsed time
elapsedTime = toc;
disp(['Elapsed time: ' num2str(elapsedTime) ' seconds']);








%========================================================================================================
% Plotting 

x1 = x(1:N+1);
x2 = x(N+2:2*N+2);
x3 = x(2*N+3:3*N+3);
[nodes,~] = LGL_nodes(N);
t0 = 0;
tf = 10;
t = ((tf-t0)/2).*nodes+(tf+t0)/2;
% t = linspace(0,10,N+1);



figure(1)
plot(t, x1);
hold on
plot(t, x2);
xlabel('Time (s)');
ylabel('State Variables');
title('Double Integrator Tracking Problem');
legend({'Positon(x1)','Velocity(x2)'},Location="northeast");
grid on
hold off


figure(2)
plot(t,x3);
xlabel('Time (s)');
ylabel('countrol variables (N)');
legend({'control variable'},Location="northeast");
title('Double integrator tracking problem');

%% Polynomial
z_value = 0.5;  % at any arbitary time in seconds
disp(['at time t =',num2str(z_value),'s']);
coeff_P = polyfit(t,x1,N);
sympref('FloatingPointOutput',true);
syms z
polynomial_p = sym(0);
for i= 1:N+1
    polynomial_p =vpa(polynomial_p+coeff_P(i)*z^(N-i+1));
end
% polynomial = polyval(coeff', z);
disp(['Position Equation:', char(polynomial_p)]);
position = subs(polynomial_p,z,z_value);
disp(['Position:', char(position),'m']);

coeff_V = polyfit(t,x2,N);
polynomial_V = sym(0);
for i= 1:N+1
    polynomial_V =vpa(polynomial_V+coeff_V(i)*z^(N-i+1));
end
disp(['Velocity Equation:', char(polynomial_V)]);
velocity = subs(polynomial_V,z,z_value);
disp(['velocity:', char(velocity),'m/s']);

coeff_A = polyfit(t,x3,N);
polynomial_A = sym(0);
for i= 1:N+1
    polynomial_A =vpa(polynomial_A+coeff_A(i)*z^(N-i+1));
end
disp(['Acceleration Equation:', char(polynomial_A)]);
accleration = subs(polynomial_A,z,z_value);
disp(['Acceleration:', char(accleration),'m/s^2']);








% Define the x and y values
x = [1, 2, 3];
y = [4, 1, 6];

% Calculate the Lagrange polynomial
lagrange_poly = lagrange_interpolation(x, y);

% Display the Lagrange polynomial
disp(lagrange_poly);

%% Polynomial
% clc;
% z_value = 0.5;
% disp(['at time t =',num2str(z_value),'s']);
% coeff_P = polyfit(t,x1,N);
% sympref('FloatingPointOutput',true);
% syms z
% polynomial_p = sym(0);
% for i= 1:N+1
%     polynomial_p =vpa(polynomial_p+coeff_P(i)*z^(N-i+1));
% end
% % polynomial = polyval(coeff', z);
% disp(['Position Equation:', char(polynomial_p)]);
% position = subs(polynomial_p,z,z_value);
% disp(['Position:', char(position),'m']);
% 
% coeff_V = polyfit(t,x2,N);
% polynomial_V = sym(0);
% for i= 1:N+1
%     polynomial_V =vpa(polynomial_V+coeff_V(i)*z^(N-i+1));
% end
% disp(['Velocity Equation:', char(polynomial_V)]);
% velocity = subs(polynomial_V,z,z_value);
% disp(['velocity:', char(velocity),'m/s']);
% 
% coeff_A = polyfit(t,x3,N);
% polynomial_A = sym(0);
% for i= 1:N+1
%     polynomial_A =vpa(polynomial_A+coeff_A(i)*z^(N-i+1));
% end
% disp(['Acceleration Equation:', char(polynomial_A)]);
% accleration = subs(polynomial_A,z,z_value);
% disp(['Acceleration:', char(accleration),'m/s^2']);

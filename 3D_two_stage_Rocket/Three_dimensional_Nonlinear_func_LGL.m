function [c,ceq,dc,dceq]= Three_dimensional_Nonlinear_func_LGL(x,M,D,problem)

dc = [];
dceq = [];


Re = problem.Re;
h_scale = problem.h_scale;
mu = problem.mu;
Omega_x = problem.Omega_x;
Omega_y = problem.Omega_y; 
Omega_z = problem.Omega_z;
rho0 = problem.rho0;
g0 = problem.g0;

m0_1 = problem.m0_1;
m0_2 = problem.m0_2;
m0 = problem.m0;
mass1_i = problem.mass1_i;
mass1_f = problem.mass1_f;
mass2_i = problem.mass2_i;
mass2_f = problem.mass2_f;

A_ref = problem.A_ref;
Cbx1 = problem.Cbx1;
Cbx2 = problem.Cbx2;
dCby1_by_dbeta = problem.dCby1_by_dbeta;
dCby2_by_dbeta = problem.dCby2_by_dbeta;
dCbz1_by_dalpha = problem.dCbz1_by_dalpha;
dCbz2_by_dalpha = problem.dCbz2_by_dalpha;

Isp = problem.Isp;
Thrust_max = problem.Thrust_max;
Thrust_max_2 = problem.Thrust_max_2;

q_max = problem.q_max;
a_sen_max = problem.a_sen_max;

hi = problem.hi;
Vi = problem.Vi;
t0 = problem.t0;

hf = problem.hf;
Vf = problem.Vf;
gamma_f = problem.gamma_f;
inclin_f = problem.inclin_f;

% Decision veriables
Rx_1 = x(1:M);
Ry_1 = x(M+1:2*M);
Rz_1 = x(2*M+1:3*M);
Vx_1 = x(3*M+1:4*M);
Vy_1 = x(4*M+1:5*M);
Vz_1 = x(5*M+1:6*M);
mass_1 = x(6*M+1:7*M);
Thrust_x1 = x(7*M+1:8*M);
Thrust_y1 = x(8*M+1:9*M);
Thrust_z1 = x(9*M+1:10*M);
Rx_2 = x(10*M+1:11*M);
Ry_2 = x(11*M+1:12*M);
Rz_2 = x(12*M+1:13*M);
Vx_2 = x(13*M+1:14*M);
Vy_2 = x(14*M+1:15*M);
Vz_2 = x(15*M+1:16*M);
mass_2 = x(16*M+1:17*M);
Thrust_x2 = x(17*M+1:18*M);
Thrust_y2 = x(18*M+1:19*M);
Thrust_z2 = x(19*M+1:20*M);
q1 = x(20*M+1:21*M);
q2 = x(21*M+1:22*M);
q3 = x(22*M+1:23*M);
q4 = x(23*M+1:24*M);
stage_time = x(24*M+1);
final_time = x(24*M+2);


% Attitude matrix
Q11 = q1.^2 - q2.^2 - q3.^2 + q4.^2;
Q12 = 2*(q1.*q2 + q3.*q4);
Q13 = 2*(q1.*q3 - q2.*q4);
Q21 = 2*(q1.*q2 - q3.*q4);
Q22 = -q1.^2 + q2.^2 - q3.^2 + q4.^2;
Q23 = 2*(q2.*q3 + q1.*q4);
Q31 = 2*(q1.*q3 + q2.*q4);
Q32 = 2*(q2.*q3 - q1.*q4);
Q33 = -q1.^2 - q2.^2 + q3.^2 + q4.^2;


Q = [Q11 Q12 Q13; Q21 Q22 Q23; Q31 Q32 Q33];

% Inertial Thrust vector
Thrust_1 = sqrt(Thrust_x1.^2 + Thrust_y1.^2 + Thrust_z1.^2);
Thrust_2 = sqrt(Thrust_x2.^2 + Thrust_y2.^2 + Thrust_z2.^2);

uTx1 = Thrust_x1./Thrust_1; 
uTy1 = Thrust_y1./Thrust_1;
uTz1 = Thrust_z1./Thrust_1;
Thrust_x1 = Thrust_1.*(Q11.*uTx1 + Q12.*uTy1 + Q13.*uTz1);
Thrust_y1 = Thrust_1.*(Q21.*uTx1 + Q22.*uTy1 + Q23.*uTz1);
Thrust_z1 = Thrust_1.*(Q31.*uTx1 + Q32.*uTy1 + Q33.*uTz1);

uTx2 = Thrust_x2./Thrust_2;
uTy2 = Thrust_y2./Thrust_2;
uTz2 = Thrust_x2./Thrust_2;
Thrust_x2 = Thrust_2.*(Q11.*uTx2 + Q12.*uTy2 + Q13.*uTz2);
Thrust_y2 = Thrust_2.*(Q21.*uTx2 + Q22.*uTy2 + Q23.*uTz2);
Thrust_z2 = Thrust_2.*(Q31.*uTx2 + Q32.*uTy2 + Q33.*uTz2);

% Gravity
g_x1 = (-mu*Rx_1)./(Rx_1.^2 + Ry_1.^2 + Rz_1.^2).^(3/2);
g_y1 = (-mu*Ry_1)./(Rx_1.^2 + Ry_1.^2 + Rz_1.^2).^(3/2);
g_z1 = (-mu*Rz_1)./(Rx_1.^2 + Ry_1.^2 + Rz_1.^2).^(3/2);
g_x2 = (-mu*Rx_2)./(Rx_2.^2 + Ry_2.^2 + Rz_2.^2).^(3/2);
g_y2 = (-mu*Ry_2)./(Rx_2.^2 + Ry_2.^2 + Rz_2.^2).^(3/2);
g_z2 = (-mu*Rz_2)./(Rx_2.^2 + Ry_2.^2 + Rz_2.^2).^(3/2);




Vrel_x1 = Vx_1 - Rz_1.*Omega_y + Ry_1.*Omega_z;
Vrel_y1 = Vy_1 - Rx_1.*Omega_z + Rz_1.*Omega_x;
Vrel_z1 = Vz_1 - Ry_1.*Omega_x + Rx_1.*Omega_y;
Vrel_x2 = Vx_2 - Rz_2.*Omega_y + Ry_2.*Omega_z;
Vrel_y2 = Vy_2 - Rx_2.*Omega_z + Rz_2.*Omega_x;
Vrel_z2 = Vz_2 - Ry_2.*Omega_x + Rx_2.*Omega_y;

Vrel_1 = sqrt(Vrel_x1.^2 + Vrel_y1.^2 + Vrel_z1.^2);
Vrel_2 = sqrt(Vrel_x2.^2 + Vrel_y2.^2 + Vrel_z2.^2);



R_1 = sqrt(Rx_1.^2 + Ry_1.^2 + Rz_1.^2);
rho_1 = rho0*exp(-(R_1-Re)./h_scale);

R_2 = sqrt(Rx_2.^2 + Ry_2.^2 + Rz_2.^2);
rho_2 = rho0*exp(-(R_2-Re)./h_scale);

Gamma = 1.4;
R = 287;
Temp0 = 288.16;
Temp_1 = Temp0 - 0.0065*(R_1-Re);
Temp_2 = Temp0 - 0.0065*(R_2-Re);

Mach_1 = Vrel_1./sqrt(Gamma*R*Temp_1);
Mach_2 = Vrel_2./sqrt(Gamma*R*Temp_2);

q_mag1 = 0.5* rho_1.* Vrel_1.^2;
q_mag2 = 0.5* rho_2.* Vrel_2.^2;


Vbx1 = Q11.*Vrel_x1 + Q21.*Vrel_y1 + Q31.*Vrel_z1;
Vby1 = Q12.*Vrel_x1 + Q22.*Vrel_y1 + Q32.*Vrel_z1;
Vbz1 = Q13.*Vrel_x1 + Q23.*Vrel_y1 + Q33.*Vrel_z1;
Vbx2 = Q11.*Vrel_x2 + Q21.*Vrel_y2 + Q31.*Vrel_z2;
Vby2 = Q12.*Vrel_x2 + Q22.*Vrel_y2 + Q32.*Vrel_z2;
Vbz2 = Q13.*Vrel_x2 + Q23.*Vrel_y2 + Q33.*Vrel_z2;



% Inertial Aerodynamic Coefficients
alpha_1 = atan(Vbz1./Vbx1);
alpha_2 = atan(Vbz2./Vbx2);
beta_1 = atan(Vby1./sqrt(Vbx1.^2 + Vbz1.^2));
beta_2 = atan(Vby2./sqrt(Vbx2.^2 + Vbz2.^2));
phi_1 = atan(Vby1./Vbx1);
phi_2 = atan(Vby2./Vbx2);

% % Cbx1 = -Ca;
% Cby1 = -Cn*sin(phi_1);
% Cbz1 = -Cn*cos(phi_1);
% % Cbx2 = -Ca;
% Cby2 = -Cn*sin(phi_2);
% Cbz2 = -Cn*cos(phi_2);

Cby1 = dCby1_by_dbeta*beta_1;
Cbz1 = dCbz1_by_dalpha*alpha_1;
Cby2 = dCby2_by_dbeta*beta_2;
Cbz2 = dCbz2_by_dalpha*alpha_2;



Cx1 = Q11.*Cbx1 + Q12.*Cby1 + Q13.*Cbz1;
Cy1 = Q21.*Cbx1 + Q22.*Cby1 + Q23.*Cbz1;
Cz1 = Q31.*Cbx1 + Q32.*Cby1 + Q33.*Cbz1;
Cx2 = Q11.*Cbx2 + Q12.*Cby2 + Q13.*Cbz2;
Cy2 = Q21.*Cbx2 + Q22.*Cby2 + Q23.*Cbz2;
Cz2 = Q31.*Cbx2 + Q32.*Cby2 + Q33.*Cbz2;


A_x1 = q_mag1.* Cx1 * A_ref;
A_y1 = q_mag1.* Cy1 * A_ref;
A_z1 = q_mag1.* Cz1 * A_ref;
A_x2 = q_mag2.* Cx2 * A_ref;
A_y2 = q_mag2.* Cy2 * A_ref;
A_z2 = q_mag2.* Cz2 * A_ref;

ceq = zeros(32*M,1);
% system dynamics
ceq(1:M,1) = D*Rx_1' - ((stage_time-t0)/2)*(Vx_1)';
ceq(M+1:2*M,1) = D*Ry_1' - ((stage_time-t0)/2)*(Vy_1)';
ceq(2*M+1:3*M,1) = D*Rz_1' - ((stage_time-t0)/2)*(Vz_1)';
ceq(3*M+1:4*M,1) = D*Vx_1' - ((stage_time-t0)/2)*((Thrust_x1 + A_x1)./mass_1 - g_x1)';
ceq(4*M+1:5*M,1) = D*Vy_1' - ((stage_time-t0)/2)*((Thrust_y1 + A_y1)./mass_1 + g_y1)';
ceq(5*M+1:6*M,1) = D*Vz_1' - ((stage_time-t0)/2)*((Thrust_z1 + A_z1)./mass_1 + g_z1)';
ceq(6*M+1:7*M,1) = D*mass_1' + ((stage_time-t0)/2)*((Thrust_1)./(g0.*Isp))';
ceq(7*M+1:8*M,1) = D*Rx_2' - ((final_time-stage_time)/2)*(Vx_2)';
ceq(8*M+1:9*M,1) = D*Ry_2' - ((final_time-stage_time)/2)*(Vy_2)';
ceq(9*M+1:10*M,1) = D*Rz_2' - ((final_time-stage_time)/2)*(Vz_2)';
ceq(10*M+1:11*M,1) = D*Vx_2' - ((final_time-stage_time)/2)*((Thrust_x2 + A_x2)./mass_2 + g_x2)';
ceq(11*M+1:12*M,1) = D*Vy_2' - ((final_time-stage_time)/2)*((Thrust_y2 + A_y2)./mass_2 + g_y2)';
ceq(12*M+1:13*M,1) = D*Vz_2' - ((final_time-stage_time)/2)*((Thrust_z2 + A_z2)./mass_2 + g_z2)';
ceq(13*M+1:14*M,1) = D*mass_2' + ((final_time-stage_time)/2)*((Thrust_2)./(g0.*Isp))';

% Initial and final constraints
ceq(14*M+1) = Rx_1(1) - (Re + hi) * cos(deg2rad(28));
ceq(14*M+2) = Vx_1(1) - Vi * cos(deg2rad(28));
ceq(14*M+3) = Ry_1(1) - 0;
ceq(14*M+4) = Vy_1(1) - 0;
ceq(14*M+5) = Rz_1(1) - (Re + hi) * sin(deg2rad(28));
ceq(14*M+6) = Vz_1(1) - Vi * sin(deg2rad(28));
ceq(14*M+7) = (Rx_1(end)^2 + Ry_1(end)^2 + Rz_1(end)^2) - ((Re + 50000)^2);
ceq(14*M+8) = (Vx_1(end)^2 + Vy_1(end)^2 + Vz_1(end)^2) - (mu/(Re + 50000));
ceq(14*M+9) = (Rx_2(end)^2 + Ry_2(end)^2 + Rz_2(end)^2) - ((Re + hf)^2);
ceq(14*M+10) = (Vx_2(end)^2 + Vy_2(end)^2 + Vz_2(end)^2) - (Vf^2);
ceq(14*M+11) = (Rx_2(end)*Vx_2(end) + Ry_2(end)*Vy_2(end) + Rz_2(end)*Vz_2(end)) - ((Re + hf) * Vf * sin(gamma_f));
ceq(14*M+12) = (Rx_2(end)*Vy_2(end) - Ry_2(end)*Vx_2(end)) - ((Re + hf) * Vf * cos(gamma_f) * sin(inclin_f));



% Knotting constraints
ceq(14*M+11) = Rx_2(1) - Rx_1(end);
ceq(14*M+12) = Ry_2(1) - Ry_1(end);
ceq(14*M+13) = Rz_2(1) - Rz_1(end);
ceq(14*M+14) = Vx_2(1) - Vx_1(end);
ceq(14*M+15) = Vy_2(1) - Vy_1(end);
ceq(14*M+16) = Vz_2(1) - Vz_1(end);
ceq(14*M+17) = mass_1(1) - mass1_i;
ceq(14*M+18) = mass_1(end) - mass1_f;
ceq(14*M+19) = mass_2(1) - mass2_i;


% Normalisation constraint for Quaternion elements 
ceq(16*M+1:17*M) = (q1.^2 + q2.^2 + q3.^2 + q4.^2) - 1;  

% % Quaternion elements
% q1 = (0.5*(1 + Q11 - Q22 - Q33).^0.5);
% ceq(17*M+1:18*M) = q2 - ((Q12 + Q21)./(4*q1));
% ceq(18*M+1:19*M) = q3 - ((Q13 + Q31)./(4*q1));
% ceq(19*M+1:20*M) = q4 - ((Q23 - Q32)./(4*q1));
% 
% ceq(20*M+1:21*M) = q2 - (0.5*(1 - Q11 + Q22 - Q33).^0.5);
% ceq(21*M+1:22*M) = q1 - ((Q21 + Q12)./(4*q2));
% ceq(22*M+1:23*M) = q3 - ((Q23 + Q32)./(4*q2));
% ceq(23*M+1:24*M) = q4 - ((Q31 - Q13)./(4*q2));
% 
% ceq(24*M+1:25*M) = q3 - (0.5*(1 - Q11 - Q22 + Q33).^0.5);
% ceq(25*M+1:26*M) = q1 - ((Q31 + Q13)./(4*q3));
% ceq(26*M+1:27*M) = q2 - ((Q32 + Q23)./(4*q3));
% ceq(27*M+1:28*M) = q4 - ((Q12 - Q21)./(4*q3));
% 
% ceq(28*M+1:29*M) = q4 - (0.5*(1 + Q11 + Q22 + Q33).^0.5);
% ceq(29*M+1:30*M) = q1 - ((Q23 - Q32)./(4*q4));
% ceq(30*M+1:31*M) = q2 - ((Q31 - Q13)./(4*q4));
% ceq(31*M+1:32*M) = q3 - ((Q12 - Q21)./(4*q4));


% Senced acceleration calculation
a_sen_x1 = D*Vx_1' - g_x1';
a_sen_y1 = D*Vy_1' - g_y1';
a_sen_z1 = D*Vz_1' - g_z1';
a_sen_mag1 = sqrt((a_sen_x1).^2 + (a_sen_y1).^2 + (a_sen_z1).^2);

a_sen_x2 = D*Vx_2' - g_x2';
a_sen_y2 = D*Vy_2' - g_y2';
a_sen_z2 = D*Vz_2' - g_z2';
a_sen_mag2 = sqrt((a_sen_x2).^2 + (a_sen_y2).^2 + (a_sen_z2).^2);


% Inequality_constraints
% c ;= []
c = zeros(10*M,1);
% % 
% c(1:M,1) = q_mag1 - q_max;
% c(M+1:2*M,1) = q_mag2 - q_max;
% 
% c(2*M+1:3*M,1) = a_sen_mag1.^2 - a_sen_max^2;
% c(3*M+1:4*M,1) = a_sen_mag2.^2 - a_sen_max^2;
% 
c(1:M,1) = Thrust_1 - Thrust_max;
c(M+1:2*M,1) = Thrust_2 - Thrust_max_2;


% Mass change constraint
for i = 1:(M-1)
    c(2*M+i) = (mass_1(i+1) - mass_1(i)) - 0;
end
for i = 1:(M-1)  
    c(2*M+(M-1)+i) = (mass_2(i+1) - mass_2(i)) - 0;
end



end


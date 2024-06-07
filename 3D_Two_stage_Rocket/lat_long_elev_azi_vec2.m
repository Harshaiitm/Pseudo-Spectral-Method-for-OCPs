function [R2_I,V2_I,qn2] = lat_long_elev_azi_vec2(M,problem) 

Re = problem.Re;
nodes = problem.nodes;
D = problem.D;

stage_time = 160;
final_time = 1680;

load time_alt_lat_long.csv
time = time_alt_lat_long(:,1);
alt = time_alt_lat_long(:,2)*1000;
lat = time_alt_lat_long(:,3);
long = time_alt_lat_long(:,4);

t_2 = ((final_time-stage_time)/2).*nodes+(final_time+stage_time)/2;

Altitude0_2 = interp1(time, alt, t_2, 'pchip');
latitude_2 = deg2rad(interp1(alt, lat, Altitude0_2, 'pchip'));
longitude_2 = deg2rad(interp1(alt, long, Altitude0_2, 'pchip'));



theta_2 = pi/2 - latitude_2;
phi_2 = longitude_2;

% Earth's rotational velocity in the inertial frame
Omega_z = 2 * pi / (24 * 60 * 60);  % Sidereal Rotation Rate (rad/s)
Omega_x = 0;
Omega_y = 0;
Omega_I = [Omega_x; Omega_y; Omega_z];


R2_B = (Altitude0_2 + Re);

% R2_B = ((R2_B(end)-R2_B(1))/2).*nodes+(R2_B(end)+R2_B(1))/2;
% theta_2 = ((theta_2(end)-theta_2(1))/2).*nodes+(theta_2(end)+theta_2(1))/2;
% phi_2 = ((phi_2(end)-phi_2(1))/2).*nodes+(phi_2(end)+phi_2(1))/2;


Rx2_I = R2_B .* sin(theta_2) .* cos(phi_2);
Ry2_I = R2_B .* sin(theta_2) .* sin(phi_2);
Rz2_I = R2_B .* cos(theta_2);

R2_I = [Rx2_I'; Ry2_I'; Rz2_I'];
R2_I_mag = sqrt(Rx2_I.^2 + Ry2_I.^2 + Rz2_I.^2);
% display(R2_I_mag);

Vx2_I = (2/(1680-160))*(D*Rx2_I);
Vy2_I = (2/(1680-160))*(D*Ry2_I);
Vz2_I = (2/(1680-160))*(D*Rz2_I);

V2_I = [Vx2_I'; Vy2_I'; Vz2_I'];
V2_I_mag = sqrt(Vx2_I.^2 + Vy2_I.^2 + Vz2_I.^2);
% disp(V2_I_mag);

Q11 = zeros(M, 1); 
Q12 = zeros(M, 1);
Q13 = zeros(M, 1);
Q21 = zeros(M, 1); 
Q22 = zeros(M, 1);
Q23 = zeros(M, 1);
Q31 = zeros(M, 1); 
Q32 = zeros(M, 1);
Q33 = zeros(M, 1);

for i = 1:M
    V_I2 = [Vx2_I(i); Vy2_I(i); Vz2_I(i)];
    V2_R = (2/(1680-160))*D(i,:)*R2_I_mag;
    V2_theta = (2/(1680-160))*R2_I_mag(i).*(D(i,:)*theta_2);
    V2_phi = (2/(1680-160))*R2_I_mag(i).*sin(theta_2(i)).*(D(i,:)*phi_2);
    V_B2 = [V2_R;V2_theta;V2_phi];    
    r = vrrotvec(V_B2, V_I2);
    Q = vrrotvec2mat(r);
    
    Q11(i) = Q(1, 1); 
    Q12(i) = Q(1, 2); 
    Q13(i) = Q(1, 3);
    Q21(i) = Q(2, 1); 
    Q22(i) = Q(2, 2); 
    Q23(i) = Q(2, 3);
    Q31(i) = Q(3, 1); 
    Q32(i) = Q(3, 2); 
    Q33(i) = Q(3, 3);
end

% Quaternion elements
q1 = -(0.5*(1 + Q11 - Q22 - Q33).^0.5);
q2 = ((Q12 + Q21)./(4*q1));
q3 = ((Q13 + Q31)./(4*q1));
q4 = ((Q23 - Q32)./(4*q1));

qn2 = [q1';q2';q3';q4'];
test = sqrt(q1.^2+q2.^2+q3.^2+q4.^2);
end






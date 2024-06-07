function [R1_I,V1_I,qn1] = lat_long_elev_azi_vec1(M,problem) 

Re = problem.Re;
nodes = problem.nodes;
D = problem.D;

t0 = 0;
stage_time = 160;

% Earth's rotational velocity in the inertial frame
Omega_z = 2 * pi / (24 * 60 * 60);  % Sidereal Rotation Rate (rad/s)
Omega_x = 0;
Omega_y = 0;
Omega_I = [Omega_x; Omega_y; Omega_z];

load time_alt_lat_long.csv
time = time_alt_lat_long(:,1);
alt = time_alt_lat_long(:,2)*1000;
lat = time_alt_lat_long(:,3);
long = time_alt_lat_long(:,4);

t_1 = ((stage_time-t0)/2).*nodes+(stage_time+t0)/2;

Altitude0_1 = interp1(time, alt, t_1, 'pchip');
latitude_1 = deg2rad(interp1(alt, lat, Altitude0_1, 'pchip'));
longitude_1 = deg2rad(interp1(alt, long, Altitude0_1, 'pchip'));

theta_1 = pi/2 - latitude_1;
phi_1 = longitude_1;

% Radial Distance from Earth's center
R1_B = (Altitude0_1 + Re);

% R1_B = ((R1_B(end)-R1_B(1))/2).*nodes+(R1_B(end)+R1_B(1))/2;
% theta_1 = ((theta_1(end)-theta_1(1))/2).*nodes+(theta_1(end)+theta_1(1))/2;
% phi_1 = ((phi_1(end)-phi_1(1))/2).*nodes+(phi_1(end)+phi_1(1))/2;

% Cartesian Coordinates
Rx1_I = R1_B .* sin(theta_1) .* cos(phi_1);
Ry1_I = R1_B .* sin(theta_1) .* sin(phi_1);
Rz1_I = R1_B .* cos(theta_1);

R1_I = [Rx1_I'; Ry1_I'; Rz1_I'];

R1_I_mag = sqrt(Rx1_I.^2 + Ry1_I.^2 + Rz1_I.^2);
% disp(R1_I_mag);

Vx1_I = (2/160)*(D*Rx1_I);
Vy1_I = (2/160)*(D*Ry1_I);
Vz1_I = (2/160)*(D*Rz1_I);

V1_I = [Vx1_I';Vy1_I'; Vz1_I'];
V1_I_mag = sqrt(Vx1_I.^2 + Vy1_I.^2 + Vz1_I.^2);


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
    V_I1 = [Vx1_I(i); Vy1_I(i); Vz1_I(i)];
    V1_R = (2/160)*D(i,:)*R1_I_mag;
    V1_theta = (2/160)*R1_I_mag(i).*(D(i,:)*theta_1);
    V1_phi = (2/160)*R1_I_mag(i).*sin(theta_1(i)).*(D(i,:)*phi_1);
    V_B1 = [V1_R;V1_theta;V1_phi];
    r = vrrotvec(V_B1, V_I1);
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
q1 = +(0.5*(1 + Q11 - Q22 - Q33).^0.5);
q2 = ((Q12 + Q21)./(4*q1));
q3 = ((Q13 + Q31)./(4*q1));
q4 = ((Q23 - Q32)./(4*q1));

qn1 = [q1';q2';q3';q4'];
test = sqrt(q1.^2+q2.^2+q3.^2+q4.^2);

end






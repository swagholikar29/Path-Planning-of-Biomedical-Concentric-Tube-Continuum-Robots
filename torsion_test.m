% syms p1 p2
clear all; clc;
dU = @tors;

x0 = [0 pi/4];

options = optimoptions('fsolve', 'OptimalityTolerance', 1e-10);
x = fsolve(dU, x0, options);
disp(rad2deg(x))


%% Webster Results

l2 = 62.3e-3;
a2 = deg2rad(265);

OD1 = 2.39e-3;
ID1 = 2.01e-3;
k1 = .0099e3;
L1 = 93.5e-3;

I1 = (pi/64)*(OD1^4 - ID1^4);
J1 = 2*I1;

OD2 = 1.6e-3;
ID2 = 0;
k2 = 0.0138e3;
L2 = 218.5e-3;

I2 = (pi/64)*(OD2^4 - ID2^4);
J2 = 2*I2;

v = .30; % .30 - .55
E = 41e6; % 41 - 75 MPa
G = E/(2*(1+v));

c1 = G*J1/L1;
c2 = G*J2/L2;
c3 = (E*I1*E*I2*k1*k2)/(E*I1+E*I2);

p1 = 1/c1*sqrt((l2*c3)^2 - ((c1*c2)/(c1+c2))^2);
p2 = p1 + acos(-c1*c2/(l2*c3*(c1+c2)));
disp(rad2deg(p1));
disp(rad2deg(p2));




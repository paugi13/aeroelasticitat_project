clc
clear
close all

c1 = 0.65;
c2 = 0.5;
h0 = 1.5;
h1 = 1.25;
h2 = 1.75;
x1 = 0.22;
x2 = 0.2;
s = 0.17;

k1 = 4000;
k2 = 3000;
k3 = 4000;
k4 = 2000;
k5 = 2500;

syms U theta1 delta theta2 theta3

% gamma1 = @(U,theta1, delta) (((-2*U*(c1-s)^2*pi*s+3*c1*U*pi*s*(c1-2))*theta1+3*c1*U*pi*s*(c1-2)*delta)/(2*(c1-s)*s+3*c1^2)+U*theta1*pi*(c1-s));
% gammaf = @(U,theta1, delta) (((-2*U*(c1-s)*pi*c1*s + 3*c1^2*U*pi*s)*theta1 + (3*c1^2*U*pi*s)*delta)/(4*(c1-s)*s+6*c1^2));
% gamma2 = @(U,theta2, theta3) ((-U*theta2 + (U*c2*theta3)/(2*(0.75*c1-0.25*c2-x1-h0+x2)))/((c2/(4*pi*(0.75*c1-0.25*c2-x1-h0+x2)*(0.75*c2-0.25*c1+x1+h0-x2)))-1/(pi*c1)));
% gamma3 = @(U,theta2, theta3) ((-U*theta3 + (U*c1*theta2)/(2*(0.75*c2-0.25*c1+x1+h0-x2)))/((c1/(4*pi*(0.75*c1-0.25*c2-x1-h0+x2)*(0.75*c2-0.25*c1+x1+h0-x2)))-1/(pi*c2)));

gamma1 =  (((-2*U*(c1-s)^2*pi*s)*theta1+3*c1*U*pi*s*(c1-s)*(theta1+delta))/(2*(c1-s)*s+3*c1^2)+U*theta1*pi*(c1-s));
gammaf =  (((-2*U*(c1-s)*pi*c1*s + 3*c1^2*U*pi*s)*theta1 + (3*c1^2*U*pi*s)*(delta))/(4*(c1-s)*s+6*c1^2));
gamma2 =  ((-U*theta2 + (U*c2*theta3)/(2*(0.75*c1-0.25*c2-x1-h0+x2)))/((c2/(4*pi*(0.75*c1-0.25*c2-x1-h0+x2)*(0.75*c2-0.25*c1+x1+h0-x2)))-1/(pi*c1)));
gamma3 =  ((-U*theta3 + (U*c1*theta2)/(2*(0.75*c2-0.25*c1+x1+h0-x2)))/((c1/(4*pi*(0.75*c1-0.25*c2-x1-h0+x2)*(0.75*c2-0.25*c1+x1+h0-x2)))-1/(pi*c2)));


C11 = 1*pi*c1;
C11_Sub1 = (3*c1*1*pi*s*(c1-s))/(2*(c1-s)*s+3*c1^2) + (3*c1^2*1*pi*s)/(4*(c1-s)*s+6*c1^2);

% C112 = (2/c1+1/s)/((4/(3*pi*c1^2))+(1/(pi*s*(c1-s))));
Cf = 3*c1*1*pi*s*(c1-s)/((2*(c1-s)*s+3*c1^2));
Cf2 = (3*c1^2*1*pi*s)/(4*(c1-s)*s+6*c1^2);
CfT = (Cf + Cf2);

C22 = -1/((c2/(4*pi*(0.75*c1-0.25*c2-x1-h0+x2)*(0.75*c2-0.25*c1+x1+h0-x2)))-1/(pi*c1));
C23 = ((1*c2)/(2*(0.75*c1-0.25*c2-x1-h0+x2)))/((c2/(4*pi*(0.75*c1-0.25*c2-x1-h0+x2)*(0.75*c2-0.25*c1+x1+h0-x2)))-1/(pi*c1));
C32 = ((1*c1)/(2*(0.75*c2-0.25*c1+x1+h0-x2)))/((c1/(4*pi*(0.75*c1-0.25*c2-x1-h0+x2)*(0.75*c2-0.25*c1+x1+h0-x2)))-1/(pi*c2));
C33 = -1/((c1/(4*pi*(0.75*c1-0.25*c2-x1-h0+x2)*(0.75*c2-0.25*c1+x1+h0-x2)))-1/(pi*c2));

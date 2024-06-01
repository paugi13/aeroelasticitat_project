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

mu = 2;

m1 = mu*h1*c1;
m2 = mu*h2*c1;
m3 = mu*h2*c2;

Is1 = (1/12)*m1*c1^2+m1*(c1/2-x1)^2;
Is2 = (1/12)*m2*c1^2+m2*(c1/2-x1)^2;
Is3 = (1/12)*m3*c2^2+m3*(c2/2-x2)^2;

% theta1 theta2 theta3 eta1 eta2
K11 = k3*x1^2;
K12 = -k3*x1^2;
K13 = 0;
K14 = 0;
K15 = 0;
K21 = -k3*x1^2;
K22 = k3*x1^2+k4*(c1-x1)^2;
K23 = 0;
K24 = (k4/h0)*(c1-x1)^2;
K25 = -(k4/h0)*(c1-x1)^2;
K31 = 0;
K32 = 0;
K33 = k5*x2^2;
K34 = k5*x2^2/h0;
K35 = -k5*x2^2/h0;
K41 = 0;
K42 = k4*(c1-x1)^2/h0;
K43 = k5*x2^2/h0;
K44 = k1 + k4*(c1-x1)^2/h0^2 + k5*x2^2/h0^2;
K45 = -k4*(c1-x1)^2/h0^2 - k5*x2^2/h0^2;
K51 = 0;
K52 = -k4*(c1-x1)^2/h0;
K53 = -k5*x2^2/h0;
K54 = -k4*(c1-x1)^2/h0^2 - k5*x2^2/h0^2;
K55 = k2 + k4*(c1-x1)^2/h0^2 + k5*x2^2/h0^2;

K = [K11 K12 K13 K14 K15;
    K21 K22 K23 K24 K25;
    K31 K32 K33 K34 K35;
    K41 K42 K43 K44 K45;
    K51 K52 K53 K54 K55];
K(2,4) = K(2,4)*-1;
K(2,5) = K(2,5)*-1;
K(3,4) = K(3,4)*-1;
K(3,5) = K(3,5)*-1;
K(4,2) = K(4,2)*-1;
K(4,3) = K(4,3)*-1;
K(5,2) = K(5,2)*-1;
K(5,3) = K(5,3)*-1;

M11 = Is1;
M12 = 0;
M13 = 0;
M14 = m1*(x1-c1/2);
M15 = 0;
M21 = 0;
M22 = Is2;
M23 = 0;
M24 = m2*(x1-c1/2);
M25 = 0;
M31 = 0;
M32 = 0;
M33 = Is3;
M34 = 0;
M35 = m3*(x2-c2/2);
M41 = m1*(x1-c1/2);
M42 = m2*(x1-c1/2);
M43 = 0;
M44 = m1 + m2;
M45 = 0;
M51 = 0;
M52 = 0;
M53 = m3*(x2-c2/2);
M54 = 0;
M55 = m3;

M = [M11 M12 M13 M14 M15;
    M21 M22 M23 M24 M25;
    M31 M32 M33 M34 M35;
    M41 M42 M43 M44 M45;
    M51 M52 M53 M54 M55];

C11 = 1*pi*c1;
C12 = 0;
C13 = 0;
C21 = 0;
C22 = -1/((c2/(4*pi*(0.75*c1-0.25*c2-x1-h0+x2)*(0.75*c2-0.25*c1+x1+h0-x2)))-1/(pi*c1));
C23 = ((1*c2)/(2*(0.75*c1-0.25*c2-x1-h0+x2)))/((c2/(4*pi*(0.75*c1-0.25*c2-x1-h0+x2)*(0.75*c2-0.25*c1+x1+h0-x2)))-1/(pi*c1));
C31 = 0;
C32 = ((1*c1)/(2*(0.75*c2-0.25*c1+x1+h0-x2)))/((c1/(4*pi*(0.75*c1-0.25*c2-x1-h0+x2)*(0.75*c2-0.25*c1+x1+h0-x2)))-1/(pi*c2));
C33 = -1/((c1/(4*pi*(0.75*c1-0.25*c2-x1-h0+x2)*(0.75*c2-0.25*c1+x1+h0-x2)))-1/(pi*c2));

U = 1;
rho = 1.25;
C = rho*[C11 C12 C13;
    C21 C22 C23; 
    C31 C32 C33];

Caux = [1 0 0 0 0;
    0 1 0 0 0;
    0 0 1 0 0;];
% l = rho*U^2*C*Caux;

S11 = h1*(x1-c1/4);
S12 = 0;
S13 = 0;
S21 = 0;
S22 = h2*(x1-c1/4);
S23 = 0;
S31 = 0;
S32 = 0;
S33 = h2*(x2-c2/4);
S41 = h1;
S42 = h2;
S43 = 0;
S51 = 0;
S52 = 0;
S53 = h2;

S = [S11 S12 S13;
    S21 S22 S23;
    S31 S32 S33;
    S41 S42 S43;
    S51 S52 S53];

A = S*C*Caux;

% APARTAT 2
neig=5;
[MODES, EIGENVAL] = eigs(K,M,neig) ;
freqs = diag(sqrt(EIGENVAL));

propertiesModes = zeros(5);

for i = 1:neig
    propertiesModes(i,1) = (MODES(4,i)+MODES(5,i))/2;
    propertiesModes(i,2) = (MODES(4,i)-MODES(5,i))/h0;
    propertiesModes(i,3) = (MODES(1,i)-MODES(2,i));
    propertiesModes(i,4) = (MODES(2,i)-propertiesModes(i,2));
    propertiesModes(i,5) = (MODES(3,i)-propertiesModes(i,2));
end

mode1 = MODES(:,1);
mode2 = MODES(:,2);
mode3 = MODES(:,3);
mode4 = MODES(:,4);
mode5 = MODES(:,5);

modeToPlot = mode5;
plot_structure(modeToPlot(4),modeToPlot(5),modeToPlot(1),modeToPlot(2),modeToPlot(3));

%APARTAT 3.
% a.
[MODES2, EIGENVAL2] = eigs(K,A,neig);
U_eig = sqrt(diag(EIGENVAL2));

% b.
Cdelta = (3*c1*1*pi*s*(c1-s))/(2*(c1-s)*s+3*c1^2) + (3*c1^2*1*pi*s)/(4*(c1-s)*s+6*c1^2);
S14 = h1*(c1-3*s/4-x1);
S24 = 0;
S34 = 0;
S44 = h1;
S54 = 0;

Sdelta = [S14;
    S24;
    S34;
    S44;
    S54];

fdelta = Cdelta*Sdelta;
D = [C11 0 0]*Caux;

U_range = 0:0.5:100; 
liftInc = zeros(length(U_range),1);

for i = 1:length(U_range)
    liftInc(i) = 1 + (U_range(i)^2/Cdelta)*D*inv(K-U_range(i)^2*A)*fdelta;
end

% figure 
% plot(U_range, liftInc, 'r')
% xlabel('U_{inf}')
% ylabel('\Delta L/L_0')
% grid on


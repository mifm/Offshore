clear all
close all

%Wind
rho = 1.225; %kg/m3 Density of air
V = 50.0; %m/s Mean Wind Speed
TI = 0.11; %At 50m/s, the turbulence intensity is 0.11. 
sigma1 = TI*V; %wind standard deviation  

%Aerodynamic properties of the blade:
Cl =  0.5; %mean lift coefficient 
Cd = 0.02; %drag coefficient throughout the blade, 
theta = 90.0; %deg relative angle
S_aero = (Cl*cos(theta) + Cd*sin(theta))^2;

f_n = 0.63;    %Hz, blade natural frequency, DTU 10MW RWT: 0.627798Hz = 1st Blade Collective Flap	
w_n = f_n*2*pi; %rad/s blade natural frequency
bl = 86.366;%m Blade length
chord  = 4; %m
chi = 0.7; %aerodynamic transfer function is 0.7. 
S_factor = rho^2 * V^2 * chord^2 * chi^2;

%Turbulence scale parameters:
Kappa1 = 42;%m, Kappa1=42m at hub height for z>42m (section 6.3 of 614100-1)
L1 = 8.1 * Kappa1; %Table B.1 of 61400-1
Lc = 8.1 * Kappa1; %coherence scale parameter 

%function of mode shape and correlation
fun = @(r1,r2)...
    (0.1351*((r1/bl).^2)+ 0.1443*((r1/bl).^3) +1.2610*((r1/bl).^4)+ 0.084*((r1/bl).^5) -0.6245*((r1/bl).^6)).*...
    (0.1351*((r2/bl).^2)+ 0.1443*((r2/bl).^3) +1.2610*((r2/bl).^4)+ 0.084*((r2/bl).^5) -0.6245*((r2/bl).^6)).*...
    exp(-12.0*((f_n.*(r1-r2)./V).^2 + (0.12*(r1-r2)./Lc).^2).^0.5);%frequency in Hz (Eq B.16 i 61400-1)

%Double interagral along the lengh of the blade:
S_I = integral2(fun,0,bl,0,bl);

%Kaimal spectra (from 614100-1 eq B.14) 
%with frequency in Hz
%NON-dimensional
S1_Kaimal_nondim = (4*sigma1^2*L1/V)/(1+6*f_n*L1/V)^(5/3);
%Dimensional S(f) is:
S1_Kaimal = S1_Kaimal_nondim*(sigma1^2)/f_n;

%???? us non-dimonsional or dimensional Kaimal?
S_FT = S_factor * S_aero * S1_Kaimal * S_I;

%Displacement
%Integrate (mass pr. length)*(mode shape^2) over the blade:
M = 1421; %kg calc done in Excel based on DTU 10MW RWT
K = (w_n^2)*M; %from presentation 6-1 slide 9+*-
%K = 2.2e9;%2.2e9Nm^2 is the flap bend stiffness at rotor radius 45m (half way between root and tip).
%Is damping ration 0.995 or 0.005 or 1.005.... ???
D = 0.005; %equal to 3% log. decr. 
Var = pi*w_n*S_FT/(K^2 * 4* D); %slide 9 of presentation 7-2. (Natural frequency (omega r) in rad/s or Hz?)
%??Why not damping squared (follows from of resuction of admittance function slide 9)
%??Where does the pi come from? 
%??is w_n rad/s or Hz?
Dev = sqrt(Var) %Tip displacement variation
%Davenprt
T=600; %sec
n_1 = sqrt(2*log(w_n*T));  %radians or Hz?
n = n_1+ 0.577/n_1;

%Extreme tip displacement
nDev = n*sqrt(Var) %m

%Plot modeshape
x = 0:0.05:1;
y = 0.1351*((x).^2)+ 0.1443*((x).^3) +1.2610*((x).^4)+ 0.084*((x).^5) -0.6245*((x).^6);

figure(1)
plot(bl*x,nDev*y,'*')
hold on
plot(bl, bl)
title('Blade mode shape')
xlabel('blade length (m)')
ylabel('mode shape (m)')
grid on









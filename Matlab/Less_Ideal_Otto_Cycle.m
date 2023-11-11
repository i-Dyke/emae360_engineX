clear
close all
clc
%{
%% Constants from Literature
eta_c = 1;                       % Compression Efficiency
eta_e = 1;                        % Expansion Efficiency
eta_cb = 0.9;                       % Combustion Efficiency
gam = 1.4;                          % Specific Heat Ratio
R_air = 0.286987;                   % Ideal Gas Const. Air (kJ/kg K)
Cp = 1.005;                         % Specific Heat for Air
Cv = R_air/(gam-1);
sub = 'AIR';                        % Substance for REFPROP calls

%% Polynomial Approx. for Cp
a1 = 3.0879/1;
a2 = (1.24597184*10^-3)/2;
a3 = (-4.23718945*10^-7)/3;
a4 = (6.74774789*10^-11)/4;
a5 = (-3.97076972*10^-15)/5;

dHi = @(Ti, To) (a1*(Ti-To) + a2*(Ti^2-To^2) + a3*(Ti^3-To^3) + a4*(Ti^4-To^4) + a5*(Ti^5-To^5))*R_air;
dHA_c = @(Ti, To) dHi(Ti, To)/eta_c;
b1_c = @(Ti, To) -(dHA_c(Ti, To)/R_air + a5*To^5 + a4*To^4 + a3*To^3 + a2*To^2 + a1*To);
TA_c = @(Ti, To) roots([a5 a4 a3 a2 a1 b1_c(Ti, To)]);

dHA_e = @(Ti, To) dHi(Ti, To)*eta_e;
b1_e = @(Ti, To) -(dHA_e(Ti, To)/R_air + a5*To^5 + a4*To^4 + a3*To^3 + a2*To^2 + a1*To);
TA_e = @(Ti, To) roots([a5 a4 a3 a2 a1 b1_e(Ti, To)]);

%% Engine Parameters
cr = 10;                             % Compression Ratio
D_tot = 1550;                       % Total Displacement, cm^3
N_cy = 4;                           % Number of Cylinders
D_cy = D_tot/N_cy;                  % Displacement per Cylinder, cm^3
rpm = 5000;                         % Engine RPM

% Conversion of volumes into meters cubed
D_tot = D_tot*10^-6;
D_cy = D_cy*10^-6;

B = 1.3;                            % Bore to Stroke Ratio
b = (4*B*D_cy/pi)^(1/3);            % Bore Size, m
s = (4*D_cy/(pi*B^2))^(1/3);        % Stroke Size, m
S_p = 2*s*rpm/60;                   % Mean Piston Speed, m/s
R = 4;                              % Connecting Rod Length to Crank Radius, m
crk_rad = s/2;                      % Crank Radius, m
con_len = crk_rad*R;                % Connecting Rod, m
MPS = 2.*s.*rpm/60;                 % Mean Piston Speed, m/s

af = 15.05;                            % Air to Fuel Ratio
Q_cb = 45000;                       % Heat from Combustion, kJ/kg

V_TDC = D_cy/(cr-1);                % Chamber Volume at TDC, m^3
V_BDC = V_TDC+D_cy;                 % Chamber Volume at BDC, m^3

V_ch = @(the) V_TDC*(1+0.5*(cr-1)*(R+1-cosd(the)-sqrt(R^2 - (sind(the)).^2)));
A_ch = @(the) (pi*b^2)/2 + (pi*b/2)*(R+1-cosd(the)-sqrt(R^2 - (sind(the)).^2))*(pi*b^2)/2;

t_deg = 60/(360*rpm);               % Time to move through 1 deg of crank, s  
t_cycle = 120/rpm;                  % Time for one thermo cycle, s

%% States

P0 = 101;                               % Initial Pressure, kPa
T0 = 293.15;                            % Initial Temperature, K
D0 = refpropm('D','T',T0,'P',P0,sub);   % Density, kg/m^3
V0 = V_TDC;

P1 = P0;
T1 = T0;
D1 = refpropm('D','T',T1,'P',P1,sub);
H1 = refpropm('H','T',T1,'P',P1,sub);   % Enthalpy, J/kg
S1 = refpropm('S','T',T1,'P',P1,sub);   % Entropy, J/kg*K
V1 = V_BDC;
m_air = V1*D1;                          % Initial mass of air, kg
fprintf("State 1:\n    Pressure = %.3f kPA\n    Temperature = %.3f K\n\n", P1, T1)

D2 = D1 * cr;
P2i = P1 * (D2/D1)^gam;
S2s = S1;
T2s = refpropm('T','P',P2i,'S',S2s,sub);
T2i = P2i/(D2*R_air);
T2 = TA_c(T2i, T1);
T2 = T2(end);
P2 = R_air*D2*T2;
H2 = refpropm('H','T',T2i,'P',P2i,sub);   % Enthalpy, J/kg
S2 = refpropm('S','T',T2i,'P',P2i,sub);   % Entropy, J/kg*K
V2 = V_TDC;
fprintf("State 2:\n    Pressure = %.3f kPA\n    Temperature = %.3f K\n\n", P2, T2)


D3 = D2;        % Combustion at constant density

% out = CEA('prob', 'UV', 'rho,kg/m^3', D2, 't,k', T2, 'o/f', af, 'fuel', 'C8H18(L)', 'C', 8, 'H', 18), 

prompt = "Enter 'i' for ideal combustion analysis, \nor the CEA combustion temperature(K): ";
% ret = input(prompt, 's');
ret = "i";
if ret == "i"
    T3 = T2i + eta_cb*(Q_cb/(Cp*af));
    P3 = D3*T3*R_air;
else
    T3 = str2double(ret);
    P3 = double(input("Enter the CEA combustion pressure(kPA): "));
end

% H3 = refpropm('H','T',T3,'P',P3,sub);   % Enthalpy, J/kg
% S3 = refpropm('S','T',T3,'P',P3,sub);   % Entropy, J/kg*K
V3 = V_TDC;
fprintf("\nState 3:\n    Pressure = %.3f kPA\n    Temperature = %.3f K\n\n", P3, T3)

D4 = D1;
P4i = P3 * (D4/D3)^gam;
T4i = P4i/(D4*R_air);
T4 = TA_e(T4i, T3);
T4 = T4(end);
P4 = R_air*D4*T4;
fprintf("State 4:\n    Pressure = %.3f kPA\n    Temperature = %.3f K\n\n", P4, T4)
V4 = V_BDC;

H4 = refpropm('H','T',T4i,'P',P4i,sub);   % Enthalpy, J/kg
S4 = refpropm('S','T',T4i,'P',P4i,sub);   % Entropy, J/kg*K

%% Outputs

W_cy = (Cv*(T3-T2i) - Cv*(T4i-T1))*m_air;         % Work Per Cylinder, kJ
W_tot = N_cy*W_cy;                              % Total Work, kJ
Po_tot = W_tot*rpm/(120);                       % Total Power, kW
Po_tot_hp = Po_tot*1.341;                       % Total Power, hp

SFC = (N_cy*m_air/af)/W_tot;                    % Specific Fuel Consumption, kg/kJ

%% Mean Piston Speed vs. Engine Speed (@ various bore-to-strokes)

figure(1)
hold on
a_B = 0.8:0.1:1.4;                          % Bore-to-strokes
a_rpm = 600:100:11000;                      % Engine Speeds, RPM

for i_B = a_B
    i_S = (4.*D_cy./(pi.*i_B))^(1/3);       % Stroke Length, m
    a_MPS = 2.*i_S.*a_rpm/60;
    plot(a_rpm, a_MPS, 'LineWidth', 1.2)
end

yline(20, "--r", 'LineWidth', 1.4)
xline(8000, "--r", 'LineWidth', 1.4)
legend([string(a_B), 'Redline MPS'],'Location', 'southeast')
title("Mean Piston Speed vs. Engine Speed (@ various bore-to-strokes)")
xlabel("Engine Speed (RPM)")
ylabel("Mean Piston Speed (m/s)")
xlim([600 11000])
ylim([0 30])

%% Compression Ratio Trade Studies

a_cr = 9:0.01:10;

figure(2);
subplot(3, 1, 1)
a_D2 = D1 .* a_cr; 
a_P2 = P1 .* (a_D2./D1).^gam;
a_T2 = a_P2./(a_D2.*R_air);

a_D3 = a_D2;
a_T3 = a_T2 + eta_cb.*(Q_cb/(Cp.*af));
a_P3 = a_D3.*a_T3.*R_air;

a_P4 = a_P3 .* (D4./a_D3).^gam;
a_T4 = a_P4./(D4.*R_air);

a_W_cy = (Cv.*(a_T3-a_T2) - Cv.*(a_T4-T1)).*m_air;
a_W_tot = N_cy.*a_W_cy;
a_Po_tot = a_W_tot.*rpm./(120) * 1.341;                       % Total Power, hp
plot(a_cr, a_Po_tot, "LineWidth", 2)
title("Total Work vs. Compression Ratio")
xlabel("Compression Ratio")
ylabel("Total Work (hp)")
ax = gca; 
ax.FontSize = 14; 


subplot(3, 1, 2)
a_SFC = (N_cy.*m_air./af)./a_W_tot;
plot(a_cr, a_SFC, "LineWidth", 2)
title("SFC vs. Compression Ratio")
xlabel("Compression Ratio")
ylabel("SFC (kJ/kg)")
ax = gca; 
ax.FontSize = 14; 

subplot(3, 1, 3)
a_fPi_tot = m_air/af * rpm/120 * 44000 *N_cy;
eff = a_Po_tot./a_fPi_tot * 100;
plot(a_cr, eff, "LineWidth", 2)
title("Efficiency vs. Compression Ratio")
xlabel("Compression Ratio")
ylabel("Efficiency (%)")
ax = gca; 
ax.FontSize = 14; 

%% Plots

figure(3)
hold on
Msize = 20;
xline(V_BDC, "--k", "LineWidth", 1.5)
xline(V_TDC, "--k", "LineWidth", 1.5)
plot(V1, P1, '.b', 'MarkerSize', Msize)
plot(V2, P2, '.b', 'MarkerSize', Msize)
plot(V3, P3, '.b', 'MarkerSize', Msize)
plot(V4, P4, '.b', 'MarkerSize', Msize)
plot(V2, P1, '.b', 'MarkerSize', Msize)
text([V2 V3 V2], [P2 P3 P1], ["State 2 > " "State 3 > " "State 0 > "], 'HorizontalAlignment', 'right')
text([V1 V4], [P1 P4], [" < State 1" " < State 4"])
ylim([-400 12500])
xlim([-0.3e-4 5e-4])
title("Pressure vs. Volume Diagram (Ideal Otto Cycle)")
xlabel("Volume (m^3)")
ylabel("Pressure (kPa)")
ax = gca; 
ax.FontSize = 10; 

%% Piston Speed

V_pis = @(the) -crk_rad .* sind(the) - ((crk_rad.^2).*sind(the).*cosd(the))/sqrt(con_len.^2 - (crk_rad*sind(the)).^2).*(2*pi*rpm/60);
the = 0:0.1:360;
plot(the, V_pis(the), 'LineWidth', 1.5)
xlim([0 360])
title({"Piston Speed Versus Crank Angle @ " + rpm + " RPM"}, {"Mean Piston Speed: " + MPS})
xlabel("Crank Angle (degrees)")
ylabel("Piston Speed (m/s)")
V_pis_max = max(V_pis(the))
yline(V_pis_max, "--r", 'LineWidth', 1.5)
yline(-V_pis_max, "--r", 'LineWidth', 1.5)
legend("Piston Speed", "Maximum Piston Speed = " + V_pis_max + " m/s", 'Location', 'southeast')
%}

%% Redefining Variables
% Basic Inputs
    % Mechanical
cr = 10;                             % Compression Ratio
D_tot = 1550;                       % Total Displacement, cm^3
N_cy = 4;                           % Number of Cylinders
D_cy = D_tot/N_cy;                  % Displacement per Cylinder, cm^3
rpm = 8000;                         % Engine Redline RPM
    % Thermo
P0 = 101;                           % Initial Pressure, kPa
T0 = 293.15;                        % Initial Temperature, K
        % (Ian CEA)
P4 = 374.14509;                     % Post Combustion (4) Pressure, kPa
T4 = 1253;                          % Post-Combustion (4) Temp, K    
        % Constants from Literature
eta_c = 1;                       % Compression Efficiency
eta_e = 1;                        % Expansion Efficiency
eta_cb = 0.9;                       % Combustion Efficiency
gam = 1.4;                          % Specific Heat Ratio
R_air = 0.286987;                   % Ideal Gas Const. Air (kJ/kg K)
Cp = 1.005;                         % Specific Heat for Air
Cv = R_air/(gam-1);


% Conversion of units
D_tot = D_tot*10^-6;
D_cy = D_cy*10^-6;
rps = rpm / 60;

b = 0.086240; %m


%% Valve Flow
clc

the = 0:0.1:720;        % Crankshaft Angle
phi = the./2;        % Camshaft Angle

TC = 0;
BC = 180;

% Valve Timings (Crank Angles)
IVO = TC-25;
IVC = BC+60;
EVO = BC-60;
EVC = 360+TC+10;

% Inlet Angle Span
IV = (IVC - IVO)/2;
EV = (EVC - EVO)/2;

L_max = 0.12*b;         % Max valve lift, m, Heywood p224
L_in = @(phi) L_max .* sind(phi/2);
L_ex = @(phi) L_max .* sind(phi/2);

% figure(4)
% polar(phi, L_in(phi))

w_in = 2e-3;               % Seat Width, m
beta_in = 30;              % Seat Angle, rad

w_ex = 2e-3;               % Seat Width, m
beta_ex = 30;              % Seat Angle, rad

D_v_in = 0.35*b;        % Inlet Valve Head Diameter, m, Heywood p222
D_v_ex = 0.28*b;       % Exhaust Valve Head Diameter, m, Heywood p222
D_v = [D_v_in; D_v_ex];

D_m_in = D_v_in - w_in;
D_m_ex = D_v_ex - w_in;
D_m = [D_m_in; D_m_ex];

D_s_in = 0.22 * D_v(1);
D_s_ex = 0.22 * D_v(2);
D_s = [D_s_in; D_s_ex];

D_p_in = 1.1 * D_v(1);
D_p_ex = 1.1 * D_v(2);
D_p = [D_p_in; D_p_ex];

L_in_0 = 0;
L_in_1 = w_in/(sind(beta_in)*cosd(beta_in));
L_in_2 = sqrt(((D_p(1)^2 - D_s(1)^2)/(4*D_m(1)))^2 - w_in^2)+w_in*tand(beta_in);
L_in_3 = L_max;

L_ex_0 = 0;
L_ex_1 = w_ex/(sind(beta_ex)*cosd(beta_ex));
L_ex_2 = sqrt(((D_p(2)^2 - D_s(2)^2)/(4*D_m(2)))^2 - w_ex^2)+w_ex*tand(beta_ex);
L_ex_3 = L_max;

A_in_1 = @(phi) pi.*L_in(phi).*cosd(beta_in).*(D_v(1)-2.*w_in+(L_in(phi)./2).*sind(2.*beta_in));
A_in_2 = @(phi) pi.*D_m(1).*sqrt((L_in(phi)-w_in.*tand(beta_in)).^2 + w_in^2);
A_in_3 = @(phi) pi/4 * (D_p(1)^2 - D_s(1)^2);

A_ex_1 = @(phi) pi.*L_ex(phi).*cosd(beta_ex).*(D_v(2)-2.*w_ex+(L_ex(phi)./2).*sind(2.*beta_ex));
A_ex_2 = @(phi) pi.*D_m(2).*sqrt((L_ex(phi)-w_ex.*tand(beta_ex)).^2 + w_ex^2);
A_ex_3 = @(phi) pi/4 * (D_p(2)^2 - D_s(2)^2);

a_L_in = L_in(phi);
[~, i_phi_in_1] = find(islocalmin(abs(a_L_in-L_in_1)));
[~, i_phi_in_2] = find(islocalmin(abs(a_L_in-L_in_2)));

a_L_ex = L_ex(phi);
[~, i_phi_ex_1] = find(islocalmin(abs(a_L_ex-L_ex_1)));
[~, i_phi_ex_2] = find(islocalmin(abs(a_L_ex-L_ex_2)));

phi_in_1 = phi(i_phi_in_1);
phi_in_2 = phi(i_phi_in_2);
i_a_phi_in = [i_phi_in_1 i_phi_in_2];
i_a_phi_in = sort(i_a_phi_in);

phi_ex_1 = phi(i_phi_ex_1);
phi_ex_2 = phi(i_phi_ex_2);
i_a_phi_ex = [i_phi_ex_1 i_phi_ex_2];
i_a_phi_ex = sort(i_a_phi_ex);

% plot(phi, L_in(phi))
% yline(L_max)
% yline(L_in_1)
% yline(L_in_2)
% xline(phi_1)
% xline(phi_2)

CD_in = 0.6;
CD_ex = 0.5;

A_in = @(phi) A_in_helper(phi, A_in_1, A_in_2, A_in_3, i_a_phi_in);
A_ex = @(phi) A_ex_helper(phi, A_ex_1, A_ex_2, A_ex_3, i_a_phi_ex);

P_in_o = P0*1.2;
P_in_T = P0;
T_in_0 = T0;    

figure(5)
in_choked = P_in_T/P_in_o <= (2/(gam+1))^(gam/(gam-1));
if in_choked
    mfr_in = @(phi) (CD_in*A_in(phi)*P_in_o)/sqrt(R_air*T_in_0) * sqrt(gam) * (2/(gam+1))^((gam+1)/(2*(gam-1)));
else
    mfr_in = @(phi) (CD_in*A_in(phi)*P_in_o)/sqrt(R_air*T_in_0) * (P_in_T/P_in_o)^(1/gam) * sqrt((2*gam)/(gam-1) * (1 - (P_in_T/P_in_o)^((gam-1)/gam)));
end
plot(phi, mfr_in(phi))
if in_choked
    title("Intake Mass Flow Rate v.s. Cam Angle (Choked)")
else
    title("Intake Mass Flow Rate v.s. Cam Angle (Unchoked)")
end
xlabel("Cam Angle (deg)")
ylabel("Mass Flow Rate (kg/s)")
xlim([0 360])
xline(phi_in_1)
xline(phi_in_2)

P_ex_o = P4;
P_ex_T = P0;
T_ex_0 = T4;

figure(6)
ex_choked = P_ex_T/P_ex_o <= (2/(gam+1))^(gam/(gam-1));
if ex_choked
    mfr_ex = @(phi) (CD_ex*A_ex(phi)*P_ex_o)/sqrt(R_air*T_ex_0) * sqrt(gam) * (2/(gam+1))^((gam+1)/(2*(gam-1)));
else
    mfr_ex = @(phi) (CD_ex*A_ex(phi)*P_ex_o)/sqrt(R_air*T_ex_0) * (P_ex_T/P_ex_o)^(1/gam) * sqrt((2*gam)/(gam-1) * (1 - (P_ex_T/P_ex_o)^((gam-1)/gam)));
end
plot(phi, mfr_ex(phi))
if in_choked
    title("Exhaust Mass Flow Rate v.s. Cam Angle (Choked)")
else
    title("Exhaust Mass Flow Rate v.s. Cam Angle (Unchoked)")
end
xlabel("Cam Angle (deg)")
ylabel("Mass Flow Rate (kg/s)")
xlim([0 360])
xline(phi_ex_1)
xline(phi_ex_2)

%% Functions

function A = A_in_helper(phi, A_in_1, A_in_2, A_in_3, i_a_phi)
    A = zeros(1, length(phi));
    A(1:i_a_phi(1)) = A_in_1(phi(1:i_a_phi(1)));
    A(i_a_phi(1):i_a_phi(2)) = A_in_2(phi(i_a_phi(1):i_a_phi(2)));
    A(i_a_phi(2):i_a_phi(3)) = A_in_3(phi(i_a_phi(2):i_a_phi(3)));
    A(i_a_phi(3):i_a_phi(4)) = A_in_2(phi(i_a_phi(3):i_a_phi(4)));
    A(i_a_phi(4):end) = A_in_1(phi(i_a_phi(4):end));
end

function A = A_ex_helper(phi, A_ex_1, A_ex_2, A_ex_3, i_a_phi)
    A = zeros(1, length(phi));
    A(1:i_a_phi(1)) = A_ex_1(phi(1:i_a_phi(1)));
    A(i_a_phi(1):i_a_phi(2)) = A_ex_2(phi(i_a_phi(1):i_a_phi(2)));
    A(i_a_phi(2):i_a_phi(3)) = A_ex_3(phi(i_a_phi(2):i_a_phi(3)));
    A(i_a_phi(3):i_a_phi(4)) = A_ex_2(phi(i_a_phi(3):i_a_phi(4)));
    A(i_a_phi(4):end) = A_ex_1(phi(i_a_phi(4):end));
end
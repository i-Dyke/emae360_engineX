clear
close all
clc

%% Constants from Literature
eta_c = 0.9;                       % Compression Efficiency
eta_e = 0.9;                        % Expansion Efficiency
eta_cb = 0.9;                       % Combustion Efficiency
gam = 1.4;                          % Specific Heat Ratio
R_air = 0.286987;                   % Ideal Gas Const. Air (kJ/kg K)
R_univ = 8.3145;                    % Universal Gas Const. (J/mol K)
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
H2 = refpropm('H','T',T2,'P',P2,sub);   % Enthalpy, J/kg
S2 = refpropm('S','T',T2,'P',P2,sub);   % Entropy, J/kg*K
V2 = V_TDC;
fprintf("State 2:\n    Pressure = %.3f kPA\n    Temperature = %.3f K\n\n", P2, T2)

D3 = D2;        % Combustion at constant density

% out = CEA('prob', 'UV', 'rho,kg/m^3', D2, 't,k', T2, 'o/f', af, 'fuel', 'C8H18(L)', 'C', 8, 'H', 18), 

prompt = "Enter 'i' for ideal combustion analysis, \nor the CEA combustion temperature(K): ";
% ret = input(prompt, 's');
ret = "i";
if ret == "i"
    T3 = T2 + eta_cb*(Q_cb/(Cp*af));
    P3 = D3*T3*R_air;
else
    T3 = 2734.1;    % K, CEA
    P3 = 9398.1;    % kPa, CEA
end

T3i = T2 + eta_cb*(Q_cb/(Cp*af));
H3 = 594.18*10^3;        % J, CEA
% H3 = refpropm('H','T',T3,'P',P3,sub);   % Enthalpy, J/kg
% S3 = refpropm('S','T',T3,'P',P3,sub);   % Entropy, J/kg*K
V3 = V_TDC;
fprintf("State 3:\n    Pressure = %.3f kPA\n    Temperature = %.3f K\n\n", P3, T3)

D4 = D1;
P4i = P3 * (D4/D3)^gam; 
T4i = P4i/(D4*R_air);
T4 = TA_e(T4i, T3);
T4 = T4(end);
P4 = R_air*D4*T4;
fprintf("State 4:\n    Pressure = %.3f kPA\n    Temperature = %.3f K\n\n", P4, T4)
V4 = V_BDC;

H4 = refpropm('H','T',T4,'P',P4,sub);   % Enthalpy, J/kg
S4 = refpropm('S','T',T4,'P',P4,sub);   % Entropy, J/kg*K

%% Outputs

W_cy = (Cv*(T3i-T2i) - Cv*(T4i-T1))*m_air;         % Work Per Cylinder, kJ
W_tot = N_cy*W_cy;                              % Total Work, kJ
Po_tot = W_tot*rpm/(120);                       % Total Power, kW
Po_tot_hp = Po_tot*1.341;                       % Total Power, hp

SFC = (N_cy*m_air/af)/W_tot;                    % Specific Fuel Consumption, kg/kJ

%% Mean Piston Speed vs. Engine Speed (@ various bore-to-strokes)

% figure(1)
% hold on
% a_B = 0.8:0.1:1.4;                          % Bore-to-strokes
% a_rpm = 600:100:11000;                      % Engine Speeds, RPM
% 
% for i_B = a_B
%     i_S = (4.*D_cy./(pi.*i_B))^(1/3);       % Stroke Length, m
%     a_MPS = 2.*i_S.*a_rpm/60;
%     plot(a_rpm, a_MPS, 'LineWidth', 1.2)
% end
% 
% yline(20, "--r", 'LineWidth', 1.4)
% xline(8000, "--r", 'LineWidth', 1.4)
% legend([string(a_B), 'Redline MPS'],'Location', 'southeast')
% title("Mean Piston Speed vs. Engine Speed (@ various bore-to-strokes)")
% xlabel("Engine Speed (RPM)")
% ylabel("Mean Piston Speed (m/s)")
% xlim([600 11000])
% ylim([0 30])

%% Compression Ratio Trade Studies

% a_cr = 9:0.01:10;
% 
% figure(2);
% subplot(3, 1, 1)
% a_D2 = D1 .* a_cr; 
% a_P2 = P1 .* (a_D2./D1).^gam;
% a_T2 = a_P2./(a_D2.*R_air);
% 
% a_D3 = a_D2;
% a_T3 = a_T2 + eta_cb.*(Q_cb/(Cp.*af));
% a_P3 = a_D3.*a_T3.*R_air;
% 
% a_P4 = a_P3 .* (D4./a_D3).^gam;
% a_T4 = a_P4./(D4.*R_air);
% 
% a_W_cy = (Cv.*(a_T3-a_T2) - Cv.*(a_T4-T1)).*m_air;
% a_W_tot = N_cy.*a_W_cy;
% a_Po_tot = a_W_tot.*rpm./(120) * 1.341;                       % Total Power, hp
% plot(a_cr, a_Po_tot, "LineWidth", 2)
% title("Total Work vs. Compression Ratio")
% xlabel("Compression Ratio")
% ylabel("Total Work (hp)")
% ax = gca; 
% ax.FontSize = 14; 
% 
% 
% subplot(3, 1, 2)
% a_SFC = (N_cy.*m_air./af)./a_W_tot;
% plot(a_cr, a_SFC, "LineWidth", 2)
% title("SFC vs. Compression Ratio")
% xlabel("Compression Ratio")
% ylabel("SFC (kJ/kg)")
% ax = gca; 
% ax.FontSize = 14; 
% 
% subplot(3, 1, 3)
% a_fPi_tot = m_air/af * rpm/120 * 44000 *N_cy;
% eff = a_Po_tot./a_fPi_tot * 100;
% plot(a_cr, eff, "LineWidth", 2)
% title("Efficiency vs. Compression Ratio")
% xlabel("Compression Ratio")
% ylabel("Efficiency (%)")
% ax = gca; 
% ax.FontSize = 14; 

%% Plots

% figure(3)
% hold on
% Msize = 20;
% xline(V_BDC, "--k", "LineWidth", 1.5)
% xline(V_TDC, "--k", "LineWidth", 1.5)
% plot(V1, P1, '.b', 'MarkerSize', Msize)
% plot(V2, P2, '.b', 'MarkerSize', Msize)
% plot(V3, P3, '.b', 'MarkerSize', Msize)
% plot(V4, P4, '.b', 'MarkerSize', Msize)
% plot(V2, P1, '.b', 'MarkerSize', Msize)
% text([V2 V3 V2], [P2 P3 P1], ["State 2 > " "State 3 > " "State 0 > "], 'HorizontalAlignment', 'right')
% text([V1 V4], [P1 P4], [" < State 1" " < State 4"])
% ylim([-400 12500])
% xlim([-0.3e-4 5e-4])
% title("Pressure vs. Volume Diagram (Ideal Otto Cycle)")
% xlabel("Volume (m^3)")
% ylabel("Pressure (kPa)")
% ax = gca; 
% ax.FontSize = 10; 

%% Piston Speed

% V_pis = @(the) -crk_rad .* sind(the) - ((crk_rad.^2).*sind(the).*cosd(the))/sqrt(con_len.^2 - (crk_rad*sind(the)).^2).*(2*pi*rpm/60);
% the = 0:0.1:360;
% plot(the, V_pis(the), 'LineWidth', 1.5)
% xlim([0 360])
% title({"Piston Speed Versus Crank Angle @ " + rpm + " RPM"}, {"Mean Piston Speed: " + MPS})
% xlabel("Crank Angle (degrees)")
% ylabel("Piston Speed (m/s)")
% V_pis_max = max(V_pis(the))
% yline(V_pis_max, "--r", 'LineWidth', 1.5)
% yline(-V_pis_max, "--r", 'LineWidth', 1.5)
% legend("Piston Speed", "Maximum Piston Speed = " + V_pis_max + " m/s", 'Location', 'southeast')

%% Valve Flow

TC = 0;
BC = 180;

% Valve Timings (Crank Angles, -360 to 360)
dwell_in = 0;
dwell_ex = 0;

IVO = TC-19;
IVC = BC+45;
EVO = BC-60;
EVC = TC+10;
IVOD = IVC-IVO;
EVOD = EVC+(360-EVO);

the_L_max_in = IVOD-dwell_in;
the_L_max_ex = EVOD-dwell_ex;

a_step = 0.01;

L_in_max = 0.12*b;          % Intake max valve lift, m, Heywood p224
L_ex_max = L_in_max*0.85;   % Exhaust max valve lift

the = a_step:a_step:the_L_max_in;
S_in = L_in_max*(64*(the./the_L_max_in).^3 - 192*(the./the_L_max_in).^4 + 192*(the./the_L_max_in).^5 - 64*(the./the_L_max_in).^6);

the = a_step:a_step:the_L_max_ex;
S_ex = L_ex_max*(64*(the./the_L_max_ex).^3 - 192*(the./the_L_max_ex).^4 + 192*(the./the_L_max_ex).^5 - 64*(the./the_L_max_ex).^6);

if dwell_in > 0
    S_in = [S_in(1:end/2), L_in_max*ones(1, dwell_in/a_step), S_in(end/2:end)];
end

if dwell_ex > 0
    S_ex = [S_ex(1:end/2), L_ex_max*ones(1, dwell_ex/a_step), S_ex(end/2:end)];
end

S_in = [zeros(1, (360+IVO)/a_step), S_in, zeros(1, (360-IVC)/a_step)]; 
S_ex = [zeros(1, EVO/a_step), S_ex, zeros(1, (360-EVC)/a_step)]; 
the = linspace(-360, 360, length(S_in));

% L_in = @(phi) L_in_helper(phi, L_in_max, phi_L_max, i_phi_L_max);
% L_ex = @(phi) L_ex_helper(phi, L_ex_max, phi_L_max, i_phi_L_max);
a_L_in = S_in;
a_L_ex = S_ex;

figure(4)
plot(the, a_L_in*1000, the, a_L_ex*1000, 'LineWidth', 2)
title("Valve Lift vs. Crank Angle")
set(gca, 'fontsize', 15);
xlabel("Crank Angle (deg)")
ylabel("Valve Lift (mm)")
xlim tight
legend("Intake", "Exhaust")

% a_L_in = [a_L_in(1:i_phi_L_max), L_in_max*ones(1, i_dwell), a_L_in(i_phi_L_max:end)];
% a_L_ex = [a_L_ex(1:i_phi_L_max), L_ex_max*ones(1, i_dwell), a_L_ex(i_phi_L_max:end)];

w_in = 2e-3;               % Seat Width, m
beta_in = 30;              % Seat Angle, rad

w_ex = 2e-3;               % Seat Width, m
beta_ex = 30;              % Seat Angle, rad

D_v_in = 0.40*b;        % Inlet Valve Head Diameter, m, Heywood p222
D_v_ex = 0.28*b;       % Exhaust Valve Head Diameter, m, Heywood p222
D_v = [D_v_in; D_v_ex];

D_m_in = D_v_in - w_in;
D_m_ex = D_v_ex - w_in;
D_m = [D_m_in; D_m_ex];

D_s_in = 0.20 * D_v(1);
D_s_ex = 0.20 * D_v(2);
D_s = [D_s_in; D_s_ex];

D_p_in = 1.1 * D_v(1);
D_p_ex = 1.1 * D_v(2);
D_p = [D_p_in; D_p_ex];

L_in_0 = 0;
L_in_1 = w_in/(sind(beta_in)*cosd(beta_in));
L_in_2 = sqrt(((D_p(1)^2 - D_s(1)^2)/(4*D_m(1)))^2 - w_in^2)+w_in*tand(beta_in);
L_in_3 = L_in_max;

L_ex_0 = 0;
L_ex_1 = w_ex/(sind(beta_ex)*cosd(beta_ex));
L_ex_2 = sqrt(((D_p(2)^2 - D_s(2)^2)/(4*D_m(2)))^2 - w_ex^2)+w_ex*tand(beta_ex);
L_ex_3 = L_ex_max;

A_in_1 = pi.*a_L_in.*cosd(beta_in).*(D_v(1)-2.*w_in+(a_L_in./2).*sind(2.*beta_in));
A_in_2 = pi.*D_m(1).*sqrt((a_L_in-w_in.*tand(beta_in)).^2 + w_in^2);
A_in_3 = pi/4 * (D_p(1)^2 - D_s(1)^2);

A_ex_1 = pi.*a_L_ex.*cosd(beta_ex).*(D_v(2)-2.*w_ex+(a_L_ex./2).*sind(2.*beta_ex));
A_ex_2 = pi.*D_m(2).*sqrt((a_L_ex-w_ex.*tand(beta_ex)).^2 + w_ex^2);
A_ex_3 = pi/4 * (D_p(2)^2 - D_s(2)^2);

[~, i_phi_in_1] = find(islocalmin(abs(a_L_in-L_in_1)));
[~, i_phi_in_2] = find(islocalmin(abs(a_L_in-L_in_2)));
if dwell_in > 0
    i_phi_in_1 = [i_phi_in_1(1) i_phi_in_1(3)];
    i_phi_in_2 = [i_phi_in_2(1) i_phi_in_2(3)];
end

[~, i_phi_ex_1] = find(islocalmin(abs(a_L_ex-L_ex_1)));
[~, i_phi_ex_2] = find(islocalmin(abs(a_L_ex-L_ex_2)));
if dwell_ex > 0
    i_phi_ex_1 = [i_phi_ex_1(1) i_phi_ex_1(3)];
    i_phi_ex_2 = [i_phi_ex_2(1) i_phi_ex_2(3)];
end

phi_in_1 = the(i_phi_in_1);
phi_in_2 = the(i_phi_in_2);
i_a_phi_in = [i_phi_in_1 i_phi_in_2];
i_a_phi_in = sort(i_a_phi_in);

phi_ex_1 = the(i_phi_ex_1);
phi_ex_2 = the(i_phi_ex_2);
i_a_phi_ex = [i_phi_ex_1 i_phi_ex_2];
i_a_phi_ex = sort(i_a_phi_ex);

% plot(phi, L_in(phi))
% yline(L_max)
% yline(L_in_1)
% yline(L_in_2)
% xline(phi_1)
% xline(phi_2)

CD_in = 0.65;
CD_ex = 0.5;

A_in = @(phi) A_in_helper(phi, A_in_1, A_in_2, A_in_3, i_a_phi_in);
A_ex = @(phi) A_ex_helper(phi, A_ex_1, A_ex_2, A_ex_3, i_a_phi_ex);

P_in_o = P0*1.1;
P_in_T = P0;
T_in_0 = T0;    

figure(5)

mfr_in_choked = @(phi) (CD_in*A_in(phi)*P_in_o)/sqrt(R_air*T_in_0) * sqrt(gam) * (2/(gam+1))^((gam+1)/(2*(gam-1)));
mfr_in = @(phi) (CD_in*A_in(phi)*P_in_o)/sqrt(R_air*T_in_0) * (P_in_T/P_in_o)^(1/gam) * sqrt((2*gam)/(gam-1) * (1 - (P_in_T/P_in_o)^((gam-1)/gam)));
a_mfr_in = mfr_in(the);
a_mfr_in_choked = mfr_in_choked(the);

plot(the, a_mfr_in, 'LineWidth', 2)
set(gca, 'fontsize', 15);
title("Intake Mass Flow Rate v.s. Crank Angle (Unchoked)")
xlabel("Crank Angle (deg)")
ylabel("Mass Flow Rate (kg/s)")
xlim tight
xline(phi_in_1(1))
xline(phi_in_1(end))
xline(phi_in_2(1)) 
xline(phi_in_2(end))

P_ex_o = P4;
P_ex_T = P0;
T_ex_0 = T4;

figure(6)
mfr_ex_choked = @(phi) (CD_ex*A_ex(phi)*P_ex_o)/sqrt(R_air*T_ex_0) * sqrt(gam) * (2/(gam+1))^((gam+1)/(2*(gam-1)));
mfr_ex = @(phi) (CD_ex*A_ex(phi)*P_ex_o)/sqrt(R_air*T_ex_0) * (P_ex_T/P_ex_o)^(1/gam) * sqrt((2*gam)/(gam-1) * (1 - (P_ex_T/P_ex_o)^((gam-1)/gam)));

a_mfr_ex = mfr_ex(the);
a_mfr_ex_choked = mfr_ex_choked(the);
plot(the, a_mfr_ex, 'LineWidth', 2)
set(gca, 'fontsize', 15);
title("Exhaust Mass Flow Rate v.s. Crank Angle (Unchoked)")

xlabel("Crank Angle (deg)")
ylabel("Mass Flow Rate (kg/s)")
xlim tight
xline(phi_ex_1(1))
xline(phi_ex_1(end))
xline(phi_ex_2(1)) 
xline(phi_ex_2(end))

% Intake SVAJ
V = diff(a_L_in*1000)/a_step * rpm * 360/60;
A = diff(V)/a_step;
J = diff(A)/a_step;

figure(8)
subplot(2, 2, 1)
sgtitle('Intake SVAJ Curves') 
plot(the, a_L_in*1000)
title("S vs. Crank Angle", "Max S = " + max(a_L_in*1000) + " mm")
xlabel("Crank Angle (deg)")
ylabel("S (mm)")
xlim tight
ylim padded

subplot(2, 2, 2)
plot(the(1:length(V)), V)
title("V vs. Crank Angle"," Max V = " + max(V) + " mm/s")
xlabel("Crank Angle (deg)")
ylabel("V (mm/s)")
xlim tight
ylim padded

subplot(2, 2, 3)
plot(the(1:length(A)), A)
title("A vs. Crank Angle", "Max A = " + max(A) + " mm/s^2")
xlabel("Crank Angle (deg)")
ylabel("A (mm/s^2)")
xlim tight
ylim padded

subplot(2,2,4)
plot(the(1:length(J)), J)
title("J vs. Crank Angle", "Max J = " + max(J) + " mm/s^3")
xlabel("Crank Angle (deg)")
ylabel("J (mm/s^3)")
xlim tight
ylim padded

% Exhuast SVAJ
V = diff(a_L_ex*1000)/a_step * rpm * 360/60;
A = diff(V)/a_step;
J = diff(A)/a_step;

figure(9)
subplot(2, 2, 1)
sgtitle('Exhaust SVAJ Curves') 
plot(the, a_L_ex*1000)
title("S vs. Crank Angle", "Max S = " + max(a_L_ex*1000) + " mm")
xlabel("Crank Angle (deg)")
ylabel("S (mm)")
xlim tight
ylim padded

subplot(2, 2, 2)
plot(the(1:length(V)), V)
title("V vs. Crank Angle"," Max V = " + max(V) + " mm/s")
xlabel("Crank Angle (deg)")
ylabel("V (mm/s)")
xlim tight
ylim padded

subplot(2, 2, 3)
plot(the(1:length(A)), A)
title("A vs. Crank Angle", "Max A = " + max(A) + " mm/s^2")
xlabel("Crank Angle (deg)")
ylabel("A (mm/s^2)")
xlim tight
ylim padded

subplot(2,2,4)
plot(the(1:length(J)), J)
title("J vs. Crank Angle", "Max J = " + max(J) + " mm/s^3")
xlabel("Crank Angle (deg)")
ylabel("J (mm/s^3)")
xlim tight
ylim padded

figure(10)
subplot(1, 2, 1)
plot(the, A_in(the)*1000000, 'LineWidth', 2)
set(gca, 'fontsize', 15);
title("Minimum Intake Flow Area vs. Crank Angle")
xlabel("Crank Angle (deg)")
ylabel("Flow Area (mm^2)")
xlim tight
xline(phi_in_1(1))
xline(phi_in_1(end))
xline(phi_in_2(1)) 
xline(phi_in_2(end))

subplot(1, 2, 2)
plot(the, A_ex(the)*1000000, 'LineWidth', 2)
set(gca, 'fontsize', 15);
title("Minimum Exhaust Flow Area vs. Crank Angle")
xlabel("Crank Angle (deg)")
ylabel("Flow Area (mm^2)")
xlim tight
xline(phi_ex_1(1))
xline(phi_ex_1(end))
xline(phi_ex_2(1)) 
xline(phi_ex_2(end))

N_in = 2;
N_ex = 2;

time = linspace(0, 2/(rpm(end)/60), length(the));
total_intake_mass = trapz(time, N_in*a_mfr_in);
min_mass = D1*V1;
a_total_mfr = N_in*a_mfr_in-N_ex*a_mfr_ex;
a_mass = zeros(1, length(time));
ins_mass = total_intake_mass+min_mass;

read_mass = 1;
file = "mass_array";
if read_mass == 0
    for i_time = 1:length(time)-1
        clc
        rng = i_time:i_time+1;
        if ins_mass < min_mass
            ins_mass = min_mass;
        else
            ins_mass = ins_mass + trapz(time(rng), a_total_mfr(rng));
        end
        a_mass(i_time) = ins_mass;
    end
    a_mass(end) = total_intake_mass+min_mass;
    ad_the = linspace(0, 720, length(the));
    % plot(ad_the, [a_mass((end/2+1):end) a_mass(1:end/2)])
    % xlim([0 720])
    % ylim([0 total_mass*1.1])
    writematrix(a_mass, file)
else
    a_mass = readmatrix(file);
end
    
%% Real Combustion
a_Vcomb = V_ch(the);

read_FHR = 1;
file = "FHR_save";
FHR = [];
if read_FHR == 0
    [FHR(1, :), FHR(2, :)] = FiniteHeatRelease(a_step, 1);
    writematrix(FHR, file)
else
    FHR = readmatrix(file);
end

a_W = FHR(2, :)*P1*V1;
a_W = a_W(1:(360/a_step));

a_P = FHR(1, :)*P1;
a_P = [a_P((end/2):end) a_P(1:(end/2-1))];
a_V = V_ch(the);
a_V(180/a_step-1) = max(a_V);
a_V(180/a_step) = max(a_V);

figure()
plot(the, FHR(1, :)*P1, 'LineWidth', 2)
set(gca, 'fontsize', 15);
title("Pressure vs. Crank Angle")
xlabel("Crank Angle (deg)")
ylabel("Pressure (kPa)")

figure()
plot(linspace(-180, 180, length(a_W)), a_W, 'LineWidth', 2)
set(gca, 'fontsize', 15);
title("Work vs. Crank Angle")
xlabel("Crank Angle (deg)")
ylabel("Work (kJ)")

figure()
plot(a_V, a_P, 'LineWidth', 2)
set(gca, 'fontsize', 15);
title("Pressure vs. Volume")
xlabel("Volume (m^3)")
ylabel("Pressure (kPa)")

% net_work = 4*FHR(2, 1)*(P1*V1);
% net_power = net_work*(rpm/120)*1.341

%% Efficiencies
close all

a_rpm = 800:8000;

a_Po_tot = zeros(1, length(a_rpm));
a_Po_tot_i = zeros(1, length(a_rpm));

a_intake_mass = zeros(1, length(a_rpm));

%a_eta_vol = -((a_rpm-3000)./20000).^2 + 0.97;
a_eta_vol = zeros(1, length(a_rpm));
a_eta_comb = zeros(1, length(a_rpm));
a_eta_thermal = zeros(1, length(a_rpm));
a_eta_mech = zeros(1, length(a_rpm)+1);
a_eta_mech(1) = 0.9;

friction_power_loss = 400*10^-3;

W_tot_real = N_cy*trapz(a_V, a_P);

max_mfr_in = max(N_in*a_mfr_in);
Q_LHV = 43.4*10^3;
T_comb = 2734.1;                                            % Peak Combustion Temp
dT = T_comb-T2;
for i_rpm = 1:length(a_rpm)
    time = linspace(0, 2/(a_rpm(i_rpm)/60), length(the));
    a_intake_mass(i_rpm) = trapz(time, N_in*a_mfr_in);
    mass = a_intake_mass(i_rpm)+min_mass;
    mass_fuel = mass/af;
    
    Q_in = mass*Cp*dT;                                      % in J
    a_eta_comb(i_rpm) = Q_in/(mass_fuel*Q_LHV);
    
    a_Po_tot_i(i_rpm) = W_tot*a_rpm(i_rpm)/(120);
    a_Po_tot(i_rpm) = W_tot_real*a_eta_vol(i_rpm)*a_eta_mech(i_rpm)*a_rpm(i_rpm)/(120);                       % Total Power, kW
    
    SFC_real = (N_cy*mass/af)/W_tot_real;                    % Specific Fuel Consumption, kg/kJ
    a_eta_thermal(i_rpm) = 1/(Q_LHV*10^3*SFC_real*10^-3);

    gross_power = a_Po_tot(i_rpm)+friction_power_loss;                         % kW
    a_eta_mech(i_rpm+1) = 1 - friction_power_loss/gross_power;
    
    % By Collin's advice, leave volumetric efficiency alone
    i_omeg = 2*pi*a_rpm(i_rpm)/60;
    %a_eta_vol(i_rpm) = 2*max_mfr_in ./ (D1*D_cy*i_omeg);
    
    if a_rpm(i_rpm) == rpm
        eta_comb = a_eta_comb(i_rpm)
        eta_thermal = a_eta_thermal(i_rpm)
        eta_mech = a_eta_mech(i_rpm+1)
        eta_vol = a_eta_vol(i_rpm)
    end
end
a_eta_mech = a_eta_mech(2:end);

% figure()
% plot(a_rpm, a_intake_mass, 'LineWidth',2)
% set(gca, 'fontsize', 13);
% title("Inducted Mass vs. Engine Speed")
% xlim padded
% ylim padded

figure()
plot(a_rpm, a_Po_tot_i, a_rpm, a_Po_tot, 'LineWidth',2)
set(gca, 'fontsize', 13);
legend("Ideal", "Real", 'Location', 'southeast')
title("Total Power vs. Engine Speed")
ylabel("Total Power (kW)")
xlabel("Engine Speed (rpm)")
xlim padded
ylim padded

figure()
plot(a_rpm, a_eta_thermal, 'LineWidth',2)
set(gca, 'fontsize', 13);
title("Thermal Efficiency vs. Engine Speed")
ylabel("Thermal Efficiency")
xlabel("Engine Speed (rpm)")
xlim padded
ylim padded

figure()
plot(a_rpm, a_eta_mech, 'LineWidth',2)
set(gca, 'fontsize', 13);
title("Mechanical Efficiency vs. Engine Speed")
ylabel("Mechanical Efficiency")
xlabel("Engine Speed (rpm)")
xlim padded
ylim([0.95 1])

figure()
plot(a_rpm, a_eta_vol, 'LineWidth',2)
set(gca, 'fontsize', 13);
title("Volumetric Efficiency (Expected) vs. Engine Speed")
ylabel("Volumetric Efficiency")
xlabel("Engine Speed (rpm)")
xlim padded
ylim([0.9 0.98])

RPM = @(R) find(a_rpm==R);

a_Po_tot(RPM(5000))
a_Po_tot(RPM(5000))*1.341

a_eta_thermal(RPM(5000))
a_eta_mech(RPM(5000))
a_eta_vol(RPM(5000))

% figure()
% a_rpm = (800+a_step):a_step:8000;
% a_eta_vol = zeros(1, length(a_rpm));
% a_A_in = A_in(the);
% CiAi = mean(a_A_in);
% A_p = pi/4 * b^2;
% a_air = 346;
% for i_rpm = 1:length(a_rpm)
%     %time = linspace(0, 2/(a_rpm(i_rpm)/60), length(the));
%     i_omeg = 2*pi*a_rpm(i_rpm)/60;
%     i_MPS = 2.*s.*a_rpm(i_rpm)/60;
%     Z = A_p*i_MPS / CiAi*a_air;
%     if Z >= 0.5
%         %i_mass = trapz(time, N_in*a_mfr_in_choked);
%         a_eta_vol(i_rpm) = 2*max_mfr_in_choked ./ (D1*D_cy*i_omeg)
%     else
%         %i_mass = trapz(time, N_in*a_mfr_in);
%         a_eta_vol(i_rpm) = 2*max_mfr_in ./ (D1*D_cy*i_omeg);
%     end
% end
% plot(a_rpm, a_eta_vol*100,'LineWidth',2)
% set(gca, 'fontsize', 15);
% title("Volumetric Efficiency vs. Engine Speed")
% xlabel("Engine Speed (RPM)")
% ylabel("Volumetric Efficiency (%)")
% xlim padded
% ylim padded

% By Collin's advice, leave volumetric efficiency alone
eta_vol = 0.9;

eta_fuel = 1/(Q_LHV*10^3*SFC*10^-3)

% Non Ideal Outputs
a1 = 3.0879;
a2 = (1.24597184*10^-3);
a3 = (-4.23718945*10^-7);
a4 = (6.74774789*10^-11);
a5 = (-3.97076972*10^-15);

%work_spc = (H3-H2)-(H4-H1);
% work_spc = 0.286987*((a5/5)*(T3^5 - T2^5)+(a4/4)*(T3^4 - T2^4)+(a3/4)*(T3^3 - T2^3)...
%             +(a2/2)*(T3^2 - T2^2)+(a1)*(T3 - T2))...
%            -0.286987*((a5/5)*(T4^5 - T1^5)+(a4/4)*(T4^5 - T1^4)+(a3/4)*(T4^3 - T1^3)...
%             +(a2/2)*(T4^2 - T1^2)+(a1)*(T4 - T1))
work_total = (Cv*(T3-T2) - Cv*(T4-T1))*(D_tot*(cr/(cr-1))*D1);
power_total = work_total*rpm*eta_mech/(2*60)
power_total_hp = power_total*1.341

%% Noise Calc

dt = 2/(rpm/60)/(720/a_step);
dP_dt = diff(a_P)/dt;
max_dP = max(dP_dt);

%% Flywheel Calculations

mass_crank = 22887.855*10^-3;
mass_piston = 491*10^-3;
mass_b = mass_piston + mass_crank/2;
sft_len = 508.5*10^-3;

omeg_min = 800*2*pi/60;
omeg_max = 8000*2*pi/60;
omeg = (omeg_max+omeg_min)/2;

a_P_torq = FHR(1, :)*P1*10^3;
a_P_torq = a_P_torq((180/a_step):(540/a_step-1));
a_P_torq = [a_P_torq a_P_torq(end)*ones(1, 360/a_step)];

shift = 180;
a_Torq_1 = a_P_torq.*crk_rad.*sind(the).*(1+(crk_rad/sft_len)*cosd(the))...
            + (mass_b./2).*(crk_rad^2).*(omeg^2)...
            .*((crk_rad/(2*sft_len))*sind(the)...
            -sind(2*the)-(3/2).*(crk_rad/sft_len).*sind(3*the));
a_Torq_2 = [a_Torq_1((end-shift/a_step):end) a_Torq_1(1:(end-shift/a_step-1))];
a_Torq_3 = [a_Torq_1((end-2*shift/a_step):end) a_Torq_1(1:(end-2*shift/a_step-1))];
a_Torq_4 = [a_Torq_1((end-3*shift/a_step):end) a_Torq_1(1:(end-3*shift/a_step-1))];
a_Torq = a_Torq_1+a_Torq_2+a_Torq_3+a_Torq_4;
mean_Torq = mean(a_Torq);
delta_E = trapz(the, a_Torq - mean_Torq);
max_rpm = 8000;
min_rpm = 800;
Cs = 2*(max_rpm-min_rpm)/(max_rpm+min_rpm);

I_flywheel = delta_E / (omeg^2*Cs);

%% Functions

function A = A_in_helper(phi, A_in_1, A_in_2, A_in_3, i_a_phi)
    A = zeros(1, length(phi));
    if length(i_a_phi) == 4
        A(1:i_a_phi(1)) = A_in_1(1:i_a_phi(1));
        A(i_a_phi(1):i_a_phi(2)) = A_in_2(i_a_phi(1):i_a_phi(2));
        A(i_a_phi(2):i_a_phi(3)) = A_in_3;
        A(i_a_phi(3):i_a_phi(4)) = A_in_2(i_a_phi(3):i_a_phi(4));
        A(i_a_phi(4):end) = A_in_1(i_a_phi(4):end);
    else
        A(1:i_a_phi(1)) = A_in_1(1:i_a_phi(1));
        A(i_a_phi(1):i_a_phi(2)) = A_in_2(i_a_phi(1):i_a_phi(2));
        A(i_a_phi(2):i_a_phi(3)) = A_in_2(i_a_phi(2):i_a_phi(3));
        A(i_a_phi(3):end) = A_in_1(i_a_phi(3):end);
    end
end

function A = A_ex_helper(phi, A_ex_1, A_ex_2, A_ex_3, i_a_phi)
    A = zeros(1, length(phi));
    if length(i_a_phi) == 4
        A(1:i_a_phi(1)) = A_ex_1(1:i_a_phi(1));
        A(i_a_phi(1):i_a_phi(2)) = A_ex_2(i_a_phi(1):i_a_phi(2));
        A(i_a_phi(2):i_a_phi(3)) = A_ex_3;
        A(i_a_phi(3):i_a_phi(4)) = A_ex_2(i_a_phi(3):i_a_phi(4));
        A(i_a_phi(4):end) = A_ex_1(i_a_phi(4):end);
    else
        A(1:i_a_phi(1)) = A_ex_1(1:i_a_phi(1));
        A(i_a_phi(1):i_a_phi(2)) = A_ex_2(i_a_phi(1):i_a_phi(2));
        A(i_a_phi(2):i_a_phi(3)) = A_ex_2(i_a_phi(2):i_a_phi(3));
        A(i_a_phi(3):end) = A_ex_1(i_a_phi(3):end);
    end
end
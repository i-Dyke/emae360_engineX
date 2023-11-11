%% Valve Flow
clc

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


the = 0:0.1:720;        % Crankshaft Angle
phi = the./2;        % Camshaft Angle

TC = 0;
BC = 180;

% Valve Timings (Crank Angles)
IVO = TC-25;
IVC = BC+60;
EVO = BC-60;
EVC = 360+TC+10;


    % CAM FUNCTION
L_max = 0.12*b;         % Max valve lift, m, Heywood p224
L_in = @(phi) L_max .* sind(phi/2);
L_ex = @(phi) L_max .* sind(phi/2);

% figure(4)
% polar(phi, L_in(phi))

    % INTAKE
w_in = 2e-3;               % Seat Width, m
beta_in = 30;              % Seat Angle, rad

    % EXHAUST
w_ex = 2e-3;               % Seat Width, m
beta_ex = 30;              % Seat Angle, rad

    % VALVE HEAD DIAMS
D_v_in = 0.35*b;        % Inlet Valve Head Diameter, m, Heywood p222
D_v_ex = 0.28*b;       % Exhaust Valve Head Diameter, m, Heywood p222
D_v = [D_v_in; D_v_ex];

    % PORT DIAMS
D_m_in = D_v_in - w_in;
D_m_ex = D_v_ex - w_in;
D_m = [D_m_in; D_m_ex];

    %
D_s_in = 0.22 * D_v(1);
D_s_ex = 0.22 * D_v(2);
D_s = [D_s_in; D_s_ex];

    %
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
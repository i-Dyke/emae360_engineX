clc
clear all
close all

format long g

D_tot = 1550;                                       % Total Displacement, cm^3
N_cy = 4;                                           % Number of Cylinders
D_cy = D_tot/N_cy;                                  % Displacement per Cylinder, cm^3
rpm = 5000;                                         % Engine RPM

% Conversion of volumes into meters cubed
D_tot = D_tot*10^-6;
D_cy = D_cy*10^-6;

cr = 10;                            % Compression Ratio
V_TDC = D_cy/(cr-1);                % Chamber Volume at TDC, m^3
V_BDC = V_TDC+D_cy;                 % Chamber Volume at BDC, m^3
V_mean = (V_TDC+V_BDC)/2;

B = 1.3;                                            % Bore to Stroke Ratio
b = (4*B*D_cy/pi)^(1/3);                            % Bore Size, m
s = (4*D_cy/(pi*B^2))^(1/3);                        % Stroke Size, m
MPS = 2*s*rpm/60;                                   % Mean Piston Speed, m/s

air = 'AIR';
T1 = 20+273.15;
T3 = 833.6069758;
T4 = 2734.1;
P1 = 101.325;
P3 = 2537.0053;
P4 = 9398.1;

k_block = 130;                                       % W/mK
k_jacket = 10;                                      % W/mK
A_block = 4*b*s*1.1;                                      % m^2
t_jacket = 3*10^-3;                                % m
t_block_1 = 15*10^-3;                                % m
t_block_2 = 18*10^-3;                               % m
R_jacket = t_jacket/(k_jacket*A_block);            % K/W
R_block_1 = t_block_1/(k_block*A_block);            % K/W
R_block_2 = t_block_2/(k_block*A_block);            % K/W

a_rpm = 800:10:8000;
a_Q_cool = zeros(1, length(a_rpm));
a_R_gas = zeros(1, length(a_rpm));

T_suf = 65+273.15;                                  % K, to prevent skin contact injury

T_cool_in = 80+273.15;                              % K, Temperature of coolant entering engine
T_cool_out = T_cool_in+25;                          % K, Temperature of coolant leaving engine (Taylor p286)
T_cool = (T_cool_in+T_cool_out)/2;

%% Combustion Chamber
eps = 0.001;                                        % Standard Surface Roughness
T_gas = (T4+T1)/2;                                  % K, average combustion temp
P_gas = (P4+P1)/2;                                  % kPa
D_gas = refpropm('D','T',T_gas,'P',P_gas,air);      % Density, kg/m^3
Pr_gas = refpropm('^','T',T_gas,'P',P_gas,air);     % Prandtl #
mu_gas = refpropm('V','T',T_gas,'P',P_gas,air);     % Dynamic Viscosity
k_gas = refpropm('L','T',T_gas,'P',P_gas,air);      % Thermal Conductivity, W/mK
L_gas = b;                                          % Characteristic Length (bore), m
mass_gas = 4.5457625e-05;                           %kg, CEA
%D_gas = mass_gas/V_TDC;

for i = 1:length(a_rpm)       
    MPS = 2*s*a_rpm(i)/60;                              % Mean Piston Speed, m/s

    V_gas = MPS;                                        % m/s
    Re_gas = D_gas*V_gas*L_gas/mu_gas;                  % Based on MPS
    
%     Nu_gas = 10.4*(Re_gas^(3/4));                       % p288 C. F. Taylor ("The Internal Combustion Engine in Theory and Practice", MIT Press, 1985)
    
    if Re_gas < 2000
        f_block = 64/Re_gas;                      % Poiseuille’s law
    else 
%       f = (0.790*log(Re_gas) - 1.64)^-2;          % Petukhov's correlation
%       f = (1.8*log10((eps/b)/3.7)^1.11 + 6.9/Re_gas)^-2;
        f_block = 0.25/(log10((eps/b)/3.7 + 5.74/Re_gas^0.9))^2;
    end
    Nu_gas = ((f_block/8)*(Re_gas-1000)*Pr_gas)/(1 + 12.7*(f_block/8)^(0.5)*(Pr_gas^(2/3) -1));
    
    h_gas = (Nu_gas*k_gas)/b;
    A_char = 4*pi/4*b^2;
    a_R_gas(i) = 1/(h_gas*A_char);

    % Air Cooling
    V_air = 20;                                         % m/s
    Q_air = air_cooling(V_air, b, s, T_suf);

    Q_gas = A_char*h_gas*(T_gas-T_cool);
    Q_cool = Q_gas - Q_air;

    a_Q_cool(i) = Q_cool;
end
figure(1)
plot(a_rpm, a_Q_cool/1000, 'LineWidth', 2)
xlim padded
ylim padded
title("Required Heat for Coolant to Disspate vs. Engine Speed")
xlabel("Engine Speed (RPM)")
ylabel("Heat to Dissipate (kW)")

%% Varying Vehicle Speed
% a_V_air = 0:0.1:35;
% a_Q_air = zeros(1, length(a_V_air));
% for j = 1:length(a_V_air)
%     a_Q_air(j) = air_cooling(a_V_air(j), b, s, k_block);
% end
% figure(2)
% plot(a_V_air, a_Q_air, 'LineWidth', 2)
% xlim padded
% ylim padded
% title("Heat Dissipated To Air vs. Vehicle Speed")
% xlabel("Vehicle Speed (m/s)")
% ylabel("Heat to Dissipate (W)")

%% Coolent Calculations
D_cool = 1055;                                      % kg/m^3, relatively insensitive to temp (https://corecheminc.com/ethylene-glycol-water-mixture-properties/#:~:text=Ethylene%20Glycol%20is%20completely%20miscible,is%20known%20as%20burst%20point.)
k_cool = 0.40;                                      % W/mK, Stone p493
mu_cool = 0.001;                                    % Dynamic Viscosity, worser case

width_cool = 15e-3;
heigh_cool = s;
A_cool = heigh_cool * width_cool;
P_cool = 2*heigh_cool + 2*width_cool;
L_cool = 4*A_cool / P_cool;

Pr_cool = 3;                                      % McCabe and Smith, 1976.
A_wall = 4*b*s;                                 % Characteristic Length (Block Length)

read = 1;
file = "coolant_velocity";
if read == 0
    cond = 0.01;
    a_V_cool = zeros(1, length(a_rpm));
    for k = 1:length(a_V_cool)
        condition = 0;
        V = 0;
        while ~condition
            Re_cool = D_cool*V*L_cool/mu_cool;

%           Nu_cool = 0.023*(Re_cool^(0.8))*(Pr_cool^(0.4));

%           syms f
%           colebrook_white = -2*log10(eps/(3.7*L_cool) + 2.51/(Re_cool*sqrt(f))) - 1/sqrt(f) == 0;
%           f = double(solve(colebrook_white, f));

            if Re_cool < 2000
                f_block = 64/Re_cool;                                 % Poiseuille’s law
%               f = (0.790*log(Re_cool) - 1.64)^-2;             % Petukhov's correlation
            else 
%               f = (100*Re_cool)^(-0.25);                           % Darcy Friction Factor, Blasius Equation
%               f = (0.790*log(Re_cool) - 1.64)^-2;             % Petukhov's correlation
%               f = (1.8*log10((eps/L_cool)/3.7)^1.11 + 6.9/Re_cool)^-2;
                f_block = 0.25/(log10((eps/L_cool)/3.7 + 5.74/Re_cool^0.9))^2;
            end
            Nu_cool = ((f_block/8)*(Re_cool-1000)*Pr_cool)/(1 + 12.7*(f_block/8)^(0.5)*(Pr_cool^(2/3) -1));

            h_cool = (Nu_cool*k_cool)/L_cool;
            R_tot = 1/(1/a_R_gas(k) + 1/R_jacket + 1/R_block_1 + 1/(h_cool*A_wall));
            T_cool_wall = T_gas - a_Q_cool(k)*R_tot;
            Q_cool = h_cool*A_wall*(T_cool_wall - T_cool);
            %(k/length(a_V_cool))*100
            if (a_Q_cool(k) - Q_cool)/a_Q_cool(k) < cond
                condition = 1;
            else
                V = V + 0.00001;
            end
        end
        a_V_cool(k) = V;
    end
    writematrix(a_V_cool, file)
else
    a_V_cool = readmatrix(file);
end

figure(3)
plot(a_rpm, a_V_cool, 'LineWidth', 2)
xlim padded
ylim padded
title("Coolant Velocity vs. Engine Speed")
xlabel("Engine Speed (rpm)")
ylabel("Coolant Velocity (m/s)")

mfr_cool = D_cool.*a_V_cool.*A_cool;
figure(4)
subplot(1, 2, 1)
plot(a_rpm, mfr_cool, 'LineWidth', 2)
xlim padded
ylim padded
title("Coolant Mass Flow Rate vs. Engine Speed")
xlabel("Engine Speed (rpm)")
ylabel("Coolant Mass Flow Rate (kg/s)")

vfr_cool = a_V_cool.*A_cool .* 1000 .* 60;
subplot(1, 2, 2)
plot(a_rpm, vfr_cool, 'LineWidth', 2)
xlim padded
ylim padded
title("Coolant Volumetric Flow Rate vs. Engine Speed")
xlabel("Engine Speed (rpm)")
ylabel("Coolant Volumetric Flow Rate (L/min)")

%% Radiator Modeling

cp_cool = 3.718*1000;                                   % J/kgK, Average cp over reasonable temperature range (https://corecheminc.com/ethylene-glycol-water-mixture-properties/#:~:text=1.084-,Specific%20Heat%C2%A0,-Specific%20Heat%2C%20or)
Q_cool_max = max(a_Q_cool);                             % W
mfr_cool_max = max(mfr_cool);
deltaT = Q_cool_max/(mfr_cool_max*cp_cool);
T_cool_out = deltaT + T_cool_in;

T_air_in = 20+273.15;                                   % K
P_air = 101.325;                                     % kPa
D_air = refpropm('D','T',T_air_in,'P',P_air,air);    % kg/m^3, Density
cp_air = refpropm('C','T',T_air_in,'P',P_air,air);   % J/kgK, Specific Heat

A_rad_ratio = 4/5;
L_rad = 400e-3;                                         % m
H_rad = L_rad*A_rad_ratio;                              % m
A_rad = L_rad*H_rad;                                    % m^2

V_air = 20;                                             % m/s

mfr_air = D_air*V_air*A_rad;                            % kg/s
T_air_out = Q_cool_max/(mfr_air*cp_air) + T_air_in;

T_air_mean = (T_air_in+T_air_out)/2;
D_air = refpropm('D','T',T_air_mean,'P',P_air,air);    % kg/m^3, Density
cp_air = refpropm('C','T',T_air_mean,'P',P_air,air);   % J/kgK, Specific Heat
Pr_air = refpropm('^','T',T_air_mean,'P',P_air,air);   % Prandtl #
k_air = refpropm('L','T',T_air_mean,'P',P_air,air);          % Thermal Conductivity, W/mK
mu_air = refpropm('V','T',T_air_mean,'P',P_air,air);         % Dynamic Viscosity

T_lm = ((T_cool_out-T_air_out)-(T_cool_in-T_air_in))/log((T_cool_out-T_air_out)/(T_cool_in-T_air_in));
P = (T_cool_out-T_cool_in)/(T_air_in-T_cool_in);
R = (T_air_in-T_air_out)/(T_cool_out-T_cool_in);
F = 0.9;                                             % Figure 11S.4, pW-42 Fund. of Heat Transfer
T_lm = F*T_lm;
UA = Q_cool_max/T_lm;                                 % W/m^2*K

% Fin Resistance Calculations
w_fin = 30e-3;
l_fin = L_rad;
t_fin = 0.1e-3;
charL_fin = w_fin + 0.5*t_fin;
P_fin = 2*l_fin + 2*t_fin;
k_fin = 237;                                           % W/mK
N_fin = floor(H_rad/(2*t_fin));
A_fin = w_fin*l_fin;
A_tot_fin = N_fin*A_fin;

Re_air = D_air*V_air*charL_fin/mu_air;
% Using flat plate correlations for Nusselt number
if Re_air < 5e5
    Nu_air = 0.664*(Re_air^(1/2))*(Pr_air^(1/3));
else
    A = 0.037*(5e5)^(4/5) - 0.664*(5e5)^(1/2);
    Nu_air = 0.037*(Re_air^(4/5) - A)*(Pr_air^(1/3));
end
h_fin = (Nu_air*k_air)/w_fin;

m_fin = sqrt((2*h_fin)/(k_fin*t_fin));
eta_fin = tanh(m_fin*charL_fin)/(m_fin*charL_fin);
eta_o = 1 - (A_fin/A_tot_fin)*(1 - eta_fin);
R_fin = 1/(eta_o*h_fin*A_tot_fin);

% Coolant Flow Resistance Calculation
Pr_cool = 7.5;                                      % McCabe and Smith, 1976.
D_pipe = 15e-3;
A_pipe = pi/4 * D_pipe^2;
V_cool = mfr_cool_max/(D_cool*A_pipe);
Re_cool = D_cool*V_cool*D_pipe/mu_cool;
if Re_cool < 2300
    Nu_cool = 4.36                         % Constant Surface Heat flux
else
    Nu_cool = 0.023*Re_cool^(4/5)*Pr_cool^(0.4);
end
h_cool = (Nu_cool*k_cool)/D_pipe;
L_pipe = 9*L_rad;
A_pipe_total = L_pipe*pi*D_pipe;
R_cool = 1/(eta_o*h_cool*A_pipe_total);

R_wall = 0;

R_tot = R_fin + R_wall + R_cool;
UA_check = 1/R_tot;

%% Pressure Drop

% Assume coolant density is constant across temperature range,
% incompressible
g = 9.81;

cool_P_res = 101.325;                       % kPa
cool_T_res = 30+273.15;                     % K

if Re_cool < 2300
    f_rad = 64/Re_cool
else
    f_rad = (1.8*log10(Re_cool/6.9))^-2;
end

D_char_rad = 4*A_pipe/(D_pipe*pi);
Ploss_rad = f_rad * (L_pipe/D_char_rad) * (V_cool^2)/2 * D_cool;

V_block = max(a_V_cool);
Re_block = D_cool*V_block*L_cool/mu_cool;

if Re_cool < 2000
    f_block = 64/Re_cool
else
    f_block = 0.25/(log10((eps/L_cool)/3.7 + 5.74/Re_cool^0.9))^2;
end

L_block = 2*4*b*1.1;
D_char_block = L_cool;
Ploss_block = f_block * (L_block/D_char_block) * (V_block^2)/2 * D_cool;

K_180 = 0.41;
K_90 = 0.50;
K_exit = 1;
K_ent = 0.5;

Ploss_minor_rad = 1/2 * D_cool * V_cool * (8*K_180 + 2*K_90);
Ploss_minor_block = 1/2 * D_cool * V_block * (4*K_90+K_exit+K_ent);

P_tot = Ploss_rad + Ploss_block + Ploss_minor_rad + Ploss_minor_block;
P_tot = P_tot/1000                     % kPa
max_vfr_cool = max(vfr_cool);
max_vfr_cool*0.26417287472922
T_cool_out-273.15
T_cool_in-273.15

close all







function Q_air = air_cooling(V_air, b, s, T_suf)
    A_block = b*s;                                      % m^2
    air = 'AIR';
    T_air = 20+273.15;                                  % K
    P_air = 101.325;                                    % kPa
    D_air = refpropm('D','T',T_air,'P',P_air,air);      % Density, kg/m^3
    Pr_air = refpropm('^','T',T_air,'P',P_air,air);     % Prandtl #
    mu_air = refpropm('V','T',T_air,'P',P_air,air);     % Dynamic Viscosity
    k_air = refpropm('L','T',T_air,'P',P_air,air);      % Thermal Conductivity, W/mK
    L_air = 4*b*1.2;                                 % Characteristic Length (Block Length)
    Re_air = D_air*V_air*L_air/mu_air;
    % Using flat plate correlations for Nusselt number
    if Re_air < 5e5
        Nu_air = 0.664*(Re_air^(1/2))*(Pr_air^(1/3));
        %Nu_air = 0.453*(Re_air^(1/2))*(Pr_air^(1/3));
    else
        Nu_air = 0.0308*(Re_air^(4/5))*(Pr_air^(1/3));
    end
    h_air = (Nu_air*k_air)/L_air;
    R_air = 1/(A_block*h_air);

    Q_air = (T_suf - T_air)/R_air;

%     w_fin = 5e-3;
%     l_fin = 1.2*(b*4);
%     t_fin = 1e-3;
%     L_fin = w_fin + 0.5*t_fin;
%     P_fin = 2*l_fin + 2*t_fin;
%     N_fin = floor(s/(3*t_fin));
% 
%     k_fin = k_block;
% 
%     A_fin = P_fin*w_fin;
%     A_tot = N_fin*A_fin + l_fin*s;
% 
%     m_fin = sqrt((2*h_air)/(k_fin*t_fin));
%     eta_fin = tanh(m_fin*L_fin)/(m_fin*L_fin);
%     eta_o = 1 - (A_tot/A_fin)*(1-eta_fin);
% 
%     Q_air = eta_o*h_air*A_tot*(T_suf - T_air);

end
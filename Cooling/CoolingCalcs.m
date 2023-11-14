clc
clear all
% close all

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
T3 = 833.6069758;
T4 = 2734.1;
P3 = 2537.0053;
P4 = 9398.1;

k_block = 130;                                       % W/mK
k_jacket = 10;                                      % W/mK
A_block = b*s;                                      % m^2
t_jacket = 3*10^-3;                                % m
t_block_1 = 12*10^-3;                                % m
t_block_2 = 18*10^-3;                               % m
R_jacket = t_jacket/(k_jacket*A_block);            % K/W
R_block_1 = t_block_1/(k_block*A_block);            % K/W
R_block_2 = t_block_2/(k_block*A_block);            % K/W

a_rpm = 800:10:8000;
a_Q_cool = zeros(1, length(a_rpm));
a_R_gas = zeros(1, length(a_rpm));

T_cool_in = 80+273.15;                              % K, Temperature of coolant entering engine
T_cool_out = T_cool_in+20;                          % K, Temperature of coolant leaving engine (Taylor p286)
T_cool = (T_cool_in+T_cool_out)/2;

%% Combustion Chamber
T_gas = (T4+T3)/2;                                  % K, average combustion temp
P_gas = (P4+P3)/2;                                  % kPa
% D_gas = refpropm('D','T',T_gas,'P',P_gas,air);      % Density, kg/m^3
Pr_gas = refpropm('^','T',T_gas,'P',P_gas,air);     % Prandtl #
mu_gas = refpropm('V','T',T_gas,'P',P_gas,air);     % Dynamic Viscosity
k_gas = refpropm('L','T',T_gas,'P',P_gas,air);      % Thermal Conductivity, W/mK
L_gas = b;                                          % Characteristic Length (bore), m
mass_gas = 4.5457625e-05;                           %kg, CEA
D_gas = mass_gas/V_mean;

for i = 1:length(a_rpm)       
    MPS = 2*s*a_rpm(i)/60;                              % Mean Piston Speed, m/s

    V_gas = MPS;                                        % m/s
    Re_gas = D_gas*V_gas*L_gas/mu_gas;                  % Based on MPS
    Nu_gas = 10.4*(Re_gas^(3/4));                         % C. F. Taylor ("The Internal Combustion Engine in Theory and Practice", MIT Press, 1985)
    h_gas = (Nu_gas*k_gas)/b;
    A_pis = pi*(b^2)/4;
    a_R_gas(i) = 1/(h_gas*A_pis);

    % Air Cooling
    V_air = 20;                                         % m/s
    Q_air = air_cooling(V_air, b, s, k_block);

    Q_gas = A_pis*h_gas*(T_gas-T_cool);
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
a_V_air = 0:0.1:35;
a_Q_air = zeros(1, length(a_V_air));
for j = 1:length(a_V_air)
    a_Q_air(j) = air_cooling(a_V_air(j), b, s, k_block);
end
figure(2)
plot(a_V_air, a_Q_air, 'LineWidth', 2)
xlim padded
ylim padded
title("Heat Dissipated To Air vs. Vehicle Speed")
xlabel("Vehicle Speed (m/s)")
ylabel("Heat to Dissipate (W)")

%% Coolent Calculations
D_cool = 1055;                                      % kg/m^3, relatively insensitive to temp (https://corecheminc.com/ethylene-glycol-water-mixture-properties/#:~:text=Ethylene%20Glycol%20is%20completely%20miscible,is%20known%20as%20burst%20point.)
k_cool = 0.47;                                      % W/mK, Stone p493
mu_cool = .0015;                                      % Dynamic Viscosity, worser case

width_cool = 15e-3;
heigh_cool = s;
A_cool = heigh_cool * width_cool;
P_cool = 2*heigh_cool + 2*width_cool;
L_cool = 4*A_cool / P_cool;

eps = 0.002;                                        % Standard Surface Roughness
Pr_cool = 135;                                      % McCabe and Smith, 1976.
A_wall = 4*b*1.1*s;                                 % Characteristic Length (Block Length)

cond = 0.01;
a_V_cool = zeros(1, length(a_rpm));
for k = 1:length(a_V_cool)
    condition = 0;
    V = 0;
    while ~condition
        clc
        Re_cool = D_cool*V*L_cool/mu_cool;
        
%       Nu_cool = 0.023*(Re_cool^(0.8))*(Pr_cool^(0.4));
        
%       syms f
%       colebrook_white = -2*log10(eps/(3.7*L_cool) + 2.51/(Re_cool*sqrt(f))) - 1/sqrt(f) == 0;
%       f = double(solve(colebrook_white, f));

        if Re_cool < 3000
            f = 64/Re_cool;                      % Poiseuille’s law
        else 
%             f = (100*Re)^(-0.25);                   % Darcy Friction Factor, Blasius Equation
            f = (0.790*log(Re_cool) - 1.64)^-2;          % Petukhov's correlation
        end
        Nu_cool = ((f/8)*(Re_cool-1000)*Pr_cool)/(1 + 12.7*(f/8)^(0.5)*(Pr_cool^(2/3) -1));
           
        h_cool = (Nu_cool*k_cool)/L_cool;
        R_tot = 1/(1/a_R_gas(k) + 1/R_jacket + 1/R_block_1 + 1/(h_cool*A_wall));
        T_cool_wall = T_gas - a_Q_cool(k)*R_tot;
        Q_cool = h_cool*A_wall*(T_cool_wall - T_cool);
        Q_cool
        a_Q_cool(k)
        abs(a_Q_cool(k) - Q_cool)/a_Q_cool(k)
        (k/length(a_V_cool))*100
        if (a_Q_cool(k) - Q_cool)/a_Q_cool(k) < cond
            condition = 1;
        else
            V = V + 0.001;
        end
    end
    V
    a_V_cool(k) = V;
end

figure(3)
plot(a_rpm, a_V_cool, 'LineWidth', 2)
xlim padded
ylim padded
title("Coolant Velocity vs. Engine Speed")
xlabel("Engine Speed (rpm)")
ylabel("Coolant Velocity (m/s)")

figure(4)
subplot(1, 2, 1)
plot(a_rpm, D_cool.*a_V_cool.*A_cool, 'LineWidth', 2)
xlim padded
ylim padded
title("Coolant Mass Flow Rate vs. Engine Speed")
xlabel("Engine Speed (rpm)")
ylabel("Coolant Mass Flow Rate (kg/s)")

subplot(1, 2, 2)
plot(a_rpm, a_V_cool.*A_cool .* 1000 .* 60, 'LineWidth', 2)
xlim padded
ylim padded
title("Coolant Volumetric Flow Rate vs. Engine Speed")
xlabel("Engine Speed (rpm)")
ylabel("Coolant Volumetric Flow Rate (L/min)")

function Q_air = air_cooling(V_air, b, s, k_block)
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
    if Re_air > 5e5
        Nu_air = 0.664*(Re_air^(1/2))*(Pr_air^(1/3));
        %Nu_air = 0.453*(Re_air^(1/2))*(Pr_air^(1/3));
    else
        Nu_air = 0.0308*(Re_air^(4/5))*(Pr_air^(1/3));
    end
    h_air = (Nu_air*k_air)/L_air;
    R_air = 1/(A_block*h_air);

    T_suf = 65+273.15;                                  % K, to prevent skin contact injury

%     Q_air_finless = (T_suf - T_air)/R_air;

    w_fin = 5e-3;
    l_fin = 1.2*(b*4);
    t_fin = 1e-3;
    L_fin = w_fin + 0.5*t_fin;
    P_fin = 2*l_fin + 2*t_fin;
    N_fin = floor(s/(3*t_fin));

    k_fin = k_block;

    A_fin = P_fin*w_fin;
    A_tot = N_fin*A_fin + l_fin*s;

    m_fin = sqrt((2*h_air)/(k_fin*t_fin));
    eta_fin = tanh(m_fin*L_fin)/(m_fin*L_fin);
    eta_o = 1 - (A_tot/A_fin)*(1-eta_fin);

    Q_air = eta_o*h_air*A_tot*(T_suf - T_air);

end
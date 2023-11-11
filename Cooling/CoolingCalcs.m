clc
clear all
close all

D_tot = 1550;                                       % Total Displacement, cm^3
N_cy = 4;                                           % Number of Cylinders
D_cy = D_tot/N_cy;                                  % Displacement per Cylinder, cm^3
rpm = 5000;                                         % Engine RPM

% Conversion of volumes into meters cubed
D_tot = D_tot*10^-6;
D_cy = D_cy*10^-6;

B = 1.3;                                            % Bore to Stroke Ratio
b = (4*B*D_cy/pi)^(1/3);                            % Bore Size, m
s = (4*D_cy/(pi*B^2))^(1/3);                        % Stroke Size, m
MPS = 2*s*rpm/60;                                   % Mean Piston Speed, m/s

air = 'AIR';
T2 = 833.6069758;z
T3 = 2734.1;

V_gas = MPS;                                        % m/s
T_gas = (T3+T2)/2;                                  % K, average combustion temp
P_gas = 9398.1;                                     % kPa
D_gas = refpropm('D','T',T_gas,'P',P_gas,air);      % Density, kg/m^3
Pr_gas = refpropm('^','T',T_gas,'P',P_gas,air);     % Prandtl #
mu_gas = refpropm('V','T',T_gas,'P',P_gas,air);     % Dynamic Viscosity
k_gas = refpropm('L','T',T_gas,'P',P_gas,air);      % Thermal Conductivity, W/mK
L_gas = b;                                          % Characteristic Length (bore)
Re_gas = D_gas*V_gas*L_gas/mu_gas;                  % Based on MPS
Nu_gas = 10.4*Re_gas^(3/4);                         % C. F. Taylor ("The Internal Combustion Engine in Theory and Practice", MIT Press, 1985)
h_gas = (Nu_gas*k_gas)/b;
A_cy = b*s;
R_gas = 1/(A_cy*h_gas);

k_block = 80;                                       % W/mK
A_block = b*s;                                      % m^2
t_block_1 = 8*10^-3;                                % m
t_block_2 = 12*10^-3;                               % m
R_block_1 = t_block_1/(k_block*A_block);            % K/W
R_block_2 = t_block_2/(k_block*A_block);            % K/W

T_air = 20+273.15;                                  % K
P_air = 101.325;                                    % kPa
D_air = refpropm('D','T',T_air,'P',P_air,air);      % Density, kg/m^3
Pr_air = refpropm('^','T',T_air,'P',P_air,air);     % Prandtl #
mu_air = refpropm('V','T',T_air,'P',P_air,air);     % Dynamic Viscosity
k_air = refpropm('L','T',T_air,'P',P_air,air);      % Thermal Conductivity, W/mK
V_air = 20;                                         % m/s
L_air = N_cy*b*1.2;                                 % Characteristic Length (Block Length)
Re_air = D_air*V_air*L_air/mu_air;
% Using flat plate correlations for Nusselt number
if Re_air > 5e5
    Nu_air = 0.0308*(Re_air^(4/5))*(Pr_air^(1/3));
else
    Nu_air = 0.453*(Re_air^(1/2))*(Pr_air^(1/3));
end
h_air = (Nu_air*k_air)/L_air;
R_air = 1/(A_block*h_air);

mfr_cool = 0.15;                                    % kg/s

D_cool = 1055;                                      % kg/m^3, relatively insensitive to temp (https://corecheminc.com/ethylene-glycol-water-mixture-properties/#:~:text=Ethylene%20Glycol%20is%20completely%20miscible,is%20known%20as%20burst%20point.)
mu_cool = 1.5;                                      % Dynamic Viscosity, worser case
A_cool = s * 6e-3;
V_cool = mfr_cool/(D_cool*A_cool);
P_cool = 2*s + 2*6e-3;
L_cool = 4*A_cool / P_cool;
Re_cool = D_cool*V_cool*L_cool/mu_cool;

R_wall = 1/(1/R_block_2 + 1/R_air);
R_cham = 1/(1/R_block_1 + 1/R_gas);

T_suf = 65+273.15;                                  % K

Q_air = (T_suf - T_air)/R_air;


clear
clc
close all
L_max=0.010348791534794;
Circle=.022;
m=265/4;
idx=(180-m)/.01;
syms S(the) V(the) A(the) J(the)
S(the) = L_max*(12.1*(the/m)^3-25.5*(the/m)^4+24.9*(the/m)^5-14.7*(the/m)^6 ...
    +4.2*(the/m)^7);
V(the) = diff(S);
A(the) = diff(V);
J(the) = diff(A);
figure(1)
the = 0.01:0.01:180;        % Crankshaft Angle
S=double(S(the));
S(idx+1:18000)=S(1:18000-idx);
S(1:idx)=0;
S=[S,flip(S)];
V=double(V(the));
V(idx+1:18000)=V(1:18000-idx);
V(1:idx)=0;
V=[V,flip(V)];
Vmax=max(V)
A=double(A(the));
A(idx+1:18000)=A(1:18000-idx);
A(1:idx)=0;
A=[A,flip(A)];
Amax=max(A)
J=double(J(the));
J(idx+1:18000)=J(1:18000-idx);
J(1:idx)=0;
J=[J,flip(J)];
Jmax=max(J)

%% 
the1 = (0.01:0.01:360);
plot(the1, S/L_max, the1, V/Vmax, the1, A/Amax, the1,J/Jmax)
ylim([-1 1])
xlim([0 360])
legend("S (m) ~ Max: " +string(round(L_max,4)), ...
    "V (m/s) ~ Max: " +string(round(Vmax*1000,4)) + "E-3", ...
    "A (m/s^2) ~ Max: " + string(round(Amax*1E6,1)) + "E-6", ...
    "J (m/s^3) ~ Max: " + string(round(Jmax*1E6,1)) + "E-6", ...
    'Location', 'southeast')


figure(2)
polarplot(the1*(pi/180),S+Circle)
title("Cam Profile 345")

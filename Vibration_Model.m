m = 114.950;             % kg, mass of block
I = 7.258;             % kg*m^2, moment of inertia of block
h = 0.35;            % m, height from CG to mount
k = 100;            % N/m, spring constant of mounting points

M = [m 0 0;...
     0 m 0;...
     0 0 I];
K = [k 0 k*h;...
     0 k 0;...
     k*h 0 k*h^2];

H = M\K;
H = sqrt(H);
f = eig(H);
omeg = 60 * sqrt(f)
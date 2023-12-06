function [a_P, a_W] = FiniteHeatRelease(steps, plots)

% Gas cycle heat release code for two engines
% engine parameters
thetas(1,1)= -10; % Engine1 start of heat release (deg)
thetas(2,1)= -10; % Engine2 start of heat release (deg)
thetad(1,1) = 40; % Engine1 duration of heat release (deg)
thetad(2,1) = 10; % Engine2 duration of heat release (deg)
r = 10; %compression ratio
gamma= 1.4; %gas const

D_cy = 1550/4 * 10^-6;
P1 = 101*10^3;                  % State 1 pressure, Pa
V1 = D_cy/(1-(1/r));            % State 1 volume, bottom dead center

m = 5.4194e-04;                 % minimum + total intake mass

R_air = 286.987;                % J/kg K
gam = 1.4;
cv = R_air/(gam-1);

T3 = 2734.1;
T2 = 736.4171;
dT = T3-T2;

Q = m*cv*dT                     % in J

q= (Q)/(P1*V1) % dimensionless total heat release Qin/P1V1

a= 5; %weibe parameter a
n= 3; %weibe exponent n
step=steps; % crankangle interval for calculation/plot
NN=360/step; % number of data points
% initialize the results data structure
save.theta=zeros(NN,1); % crankangle
save.vol=zeros(NN,1); % volume
save.press=zeros(NN,2); % pressure
save.work=zeros(NN,2); % work
pinit(1) = P1/P1; % Engine 1 initial dimensionless pressure P/P1
pinit(2) = P1/P1; % Engine 2 initial dimensionless pressure P/P1
% for loop for engine1 and engine2
for j=1:2
    theta = -180; %initial crankangle
    thetae = theta + step; %final crankangle in step
    fy(1) = pinit(j); % assign initial pressure to working vector
    fy(2) =0.; % reset work vector
    % for loop for pressure and work calculation
    for i=1:NN
        [fy, vol] = integrate(theta,thetae,fy);
        % reset to next interval
        theta = thetae;
        thetae = theta+step;
        % copy results to output vectors
        save.theta(i)=theta;
        save.vol(i)=vol;
        save.press(i,j)=fy(1);
        save.work(i,j)=fy(2);
    end %end of pressure and work iteration loop
end %end of engine iteration loop
[pmax1, id_max1] = max(save.press(:,1)); %Engine 1 max pressure
[pmax2, id_max2] = max(save.press(:,2)); %Engine 2 max pressure
thmax1=save.theta(id_max1);%Engine 1 crank angle
thmax2=save.theta(id_max2);%Engine 2 crank angle
w1=save.work(NN,1);
w2=save.work(NN,2);
eta1= w1/q; % thermal efficiency
eta2= w2/q;
imep1 = eta1*q*(r/(r -1)); %imep
imep2 = eta2*q*(r/(r -1));
eta_rat1 = eta1/(1-r^(1-gamma));
eta_rat2 = eta2/(1-r^(1-gamma));
% output overall results
fprintf(' Engine 1 Engine 2 \n');
fprintf(' Theta_start %5.2f %5.2f \n', thetas(1,1), thetas(2,1));
fprintf(' Theta_dur %5.2f %5.2f \n', thetad(1,1), thetad(2,1));
fprintf(' P_max/P_1 %5.2f %5.2f \n', pmax1, pmax2);
fprintf(' Theta_max %7.1f %7.1f \n',thmax1,thmax2);
fprintf(' Net Work/P1V1 %7.2f %7.2f \n', w1,w2);
fprintf(' Efficiency %5.3f %5.3f \n', eta1, eta2);
fprintf(' Eff. Ratio %5.3f %5.3f \n', eta_rat1, eta_rat2);
fprintf(' Imep/P1 %5.2f %5.2f \n', imep1, imep2);
%plot results
if plots
    figure()
    plot(save.theta,save.press(:,1),'-','linewidth',2 )
    set(gca, 'fontsize', 18,'linewidth',2);
%     legend('Engine 1', 'Engine 2','Location','NorthWest')
    xlabel('Theta (deg)','fontsize', 18)
    ylabel('Pressure (bar)','fontsize', 18)
    print -deps2 heatrelpressure
    figure()
    plot(save.theta,save.work(:,1),'-', 'linewidth',2)
    set(gca, 'fontsize', 18,'linewidth',2);
%     legend('Engine 1', 'Engine 2','Location','NorthWest')
    xlabel('Theta (deg)','fontsize', 18)
    ylabel('Work','fontsize', 18)
end

function[fy,vol] = integrate(theta,thetae,fy)
    %ode23 integration of the pressure differential equation
    %from theta to thetae with current values of fy as initial conditions
    [tt, yy] = ode23(@rates, [theta thetae], fy);
    for k=1:2
        fy(k) = yy(length(tt),k); %put last element of yy into fy vector
    end
    %pressure differential equation
    function [yprime] = rates(theta,fy)
        vol=(1.+ (r -1)/2.*(1-cosd(theta)))/r;
        dvol=(r - 1)/2.*sind(theta)/r*pi/180.; %dvol/dtheta
        dx=0.; %set heat release to zero
        if(theta > thetas(j)) % then heat release dx > 0
            dum1=(theta -thetas(j))/thetad(j);
            x=1.- exp(-(a*dum1^n));
            dx=(1-x)*a*n*dum1^(n-1)/thetad(j); %dx/dthetha
        end
        term1= -gamma*fy(1)*dvol/vol;
        term2= (gamma-1)*q*dx/vol;
        yprime(1,1)= term1 + term2;
        yprime(2,1)= fy(1)*dvol;
    end %end of function rates
end %end of function integrate2
a_P = save.press(:,1)';
a_P = [ones(1, (180/step)), a_P, a_P(end)*ones(1, (180/step))];

a_W = save.work(:,1)';
a_W = [a_W, zeros(1, length(a_P)-length(a_W))];
end % heat_release_weibe2
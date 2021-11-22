% Spacecraft Guidance and Navigation 2021/2022
% Assignment # 1
% Author: Niccolò Bardazzi 10800456

%% Ex 1
clearvars; close all; clc; 

attractor = nbody_init({'Earth'});
mu = attractor{1}.GM;

% Point 1, validation:
r0 = [29597.43; 0; 0];    
v0 = cross([0; 0; 1], r0)/norm(cross([0; 0; 1], r0))*sqrt(mu/norm(r0));
x0 = [r0; v0];
T_start = 2*pi*sqrt(norm(r0)^3/mu);
options_ode = odeset('reltol', 1e-13, 'abstol', 1e-16); % high tolerances just for better drawing
[~, y_start] = ode113(@(t,x) ode_orbit(t,x,mu), [0 T_start], x0, options_ode); 
fprintf('Error on integration is: %.9e\n', norm(y_start(end, 1:3)-y_start(1, 1:3)))

% Point 2
x0_guess = [r0; v0];
x_target = [-10000; 8500; 10];
data.x_target = x_target;
data.mu = mu;
options_ode = odeset('reltol', 1e-6, 'abstol', 1e-6);
tof = T_start/2;
tol = 1.e-4;
% Shooting method with STM
[STM.x0, STM.x_target, STM.count]  = NM_STM(x0_guess, x_target, mu, [T_start T_start+tof], tol, options_ode);
% Shooting with finite differences
[FD.x0, FD.x_target, FD.count] = NM_FD(x0_guess, x_target, mu, [T_start T_start+tof], tol, options_ode);
% Clasic Lambert's solver
[~, ~, ~, ~, VI, ~, ~, ~] = lambertMR(r0, x_target, tof, mu, 0, 0, 0, 0);
fprintf('Error of shooting with STM with respect to the classic Lambert''s solver is:\n %.8e \n',norm(VI(1:3)'-STM.x0(4:6)));
fprintf('Error of shooting with FD with respect to the classic Lambert''s solver is:\n %.8e \n %.8e \n %.8e \n',norm((VI(1:3)'-FD.x0(4:6))));

figure(1)
options = odeset('reltol', 1e-13, 'abstol', 1e-16, 'event', @target_reached);
[~, y, te, ye, ie] = ode113(@ode_orbit_target, [T_start T_start+tof], STM.x0, options, data); 
scatter3(x0(1), x0(2),x0(3),'LineWidth',1.3), hold on
scatter3(x_target(1), x_target(2), x_target(3),2.e2,'Marker','x','LineWidth',1.3), hold on
plot3(y(:,1), y(:,2), y(:,3), 'LineStyle', '-.', 'Color','k','LineWidth',1), hold on
v_target = cross([0; 0; 1], x_target)/norm(cross([0; 0; 1], x_target))*sqrt(mu/norm(x_target));
T_target = 2*pi*sqrt(norm(x_target)^3/mu);
[~, y_target] = ode113(@(t,x) ode_orbit(t,x,mu), [0 T_target], [x_target; v_target], options_ode);
title('Lambert''s problem')
RE = cspice_bodvrd('Earth','RADII',3);
[X, Y, Z] = sphere;
surf(X*RE(1), Y*RE(2), Z*RE(3)), hold off
axis equal
txt = {'Start'};
text(x0(1)-1e3,x0(2)+1e3,x0(3)-2.5e3,txt,'FontSize',15)
txt = {'Target'};
text(x_target(1)-1e3,x_target(2)+1e3,x_target(3)-2.5e3,txt,'FontSize',15)
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')


% Point 3
labels = {'Sun';
          'Earth';
          'Mars Barycenter';
          };

body = nbody_init(labels(2:end));
attractor = nbody_init(labels(1));
frame = 'ECLIPJ2000';
mu = attractor{1}.GM;

% Adimensionalization
LL = 1e7; %[km]
TT = 1e7; %[s]
GMGM = (LL^3)/(TT^2);

mu = mu/GMGM;
AU = cspice_convrt(1,'AU','km');

min_dep = cspice_str2et('2033-Apr-1 00:00:00.0000')./TT; 
max_dep = cspice_str2et('2036-Apr-1 00:00:00.0000')./TT;   
min_arr = cspice_str2et('2033-Sep-1 00:00:00.0000')./TT;
max_arr = cspice_str2et('2036-Sep-1 00:00:00.0000')./TT; 


tspan = linspace(min_dep, max_arr, 1000).*TT;

dep_guess_dim = cspice_spkezr('EARTH', cspice_str2et('Oct 1 00:00:00 UTC 2035'), frame, 'NONE', 'SUN');
% dep_guess_dim = cspice_spkezr('EARTH', cspice_str2et('Oct 13 00:00:00 UTC 2033'), frame, 'NONE', 'SUN');
dep_guess = [dep_guess_dim(1:3)./LL; dep_guess_dim(4:6)./(LL/TT)]; % adimensionalized

x1_guess = cspice_spkezr(body{1}.name, min_dep, frame, 'NONE', 'SSB');
x0 = [cspice_str2et('Oct 1 00:00:00 UTC 2035')/TT; 
      cspice_str2et('May 20 00:00:00 UTC 2036')/TT; 
      dep_guess(1:3)+[10;10;0]; 
      dep_guess(4:6)+[10;1;5]]; 

% x0 = [cspice_str2et('Oct 13 00:00:00 UTC 2033')/TT; 
%       cspice_str2et('Jun 20 00:00:00 UTC 2034')/TT; 
%       dep_guess(1:3)+[10;10;0]; 
%       dep_guess(4:6)+[10;1;5]]; 

% fmincon set
A = [1,-1,0,0,0,0,0,0];  b = 0;   
Aeq = []; beq = [];
lb = [min_dep, min_arr, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf];
ub = [max_dep, max_arr, Inf, Inf, Inf, Inf, Inf, Inf];
data.body{1} = body{1};   data.body{2} = body{2}; 
data.frame = frame;
data.mu = mu;
data.TT = TT;               data.LL = LL;            data.GMGM = GMGM;
options = optimset('Display', 'iter', 'Algorithm', 'interior-point', 'TolFun', 1.e-4, 'TolCon', 1.e-8, 'Maxiter', 250);
on = 1; off = 0;
data.turn = on;

if data.turn == on
    figure(1)
    xE0 = cspice_spkezr(body{1}.name, min_dep.*TT, frame, 'NONE', 'SSB');
    xM0 = cspice_spkezr(body{2}.name, min_dep.*TT, frame, 'NONE', 'SSB');
    options_ode = odeset('reltol', 1e-9, 'abstol', 1e-9);
    T_Eprop = 688*24*cspice_spd;
    [~, xEprop] = ode113(@(t,x) ode_orbit(t,x,mu.*GMGM), [0 T_Eprop], xE0, options_ode); 
    T_Mprop = 366*24*cspice_spd;
    [~, xMprop] = ode113(@(t,x) ode_orbit(t,x,mu.*GMGM), [0 T_Mprop], xM0, options_ode); 
    plot3(xEprop(:,1), xEprop(:,2), xEprop(:,3),'color',[0, 0.5, 0],'LineStyle','-','LineWidth',1), hold on
    plot3(xMprop(:,1), xMprop(:,2), xMprop(:,3),'color',[0.8500, 0.3250, 0.0980],'LineStyle','-','LineWidth',1), hold on  
    xlabel('x [km]')
    ylabel('y [km]')
    zlabel('z [km]')
    title(sprintf('Iteration process'),'Interpreter','latex','FontSize',16)
end

[y, dv] = fmincon(@obj, x0, A, b, Aeq, beq, lb, ub, @nonlcons, options, data);

if data.turn == on
    xE0 = cspice_spkezr(body{1}.name, y(1).*TT, frame, 'NONE', 'SSB');
    xM0 = cspice_spkezr(body{2}.name, y(2).*TT, frame, 'NONE', 'SSB');
    RE = cspice_bodvrd('Earth','RADII',3);
    [X, Y, Z] = sphere;
    surf(X*RE(1)+xE0(1), Y*RE(2)+xE0(2), Z*RE(3)+xE0(3)), hold on
    txt = {'Earth'};
    text(xE0(1)-1.5e7,xE0(2)-1.5e7,xE0(3)-1.e6,txt,'FontSize',14)
    RM = cspice_bodvrd('Mars','RADII',3);
    surf(X*RM(1)+xM0(1), Y*RM(2)+xM0(2), Z*RM(3)+xM0(3)), hold on
    txt = {'Mars'};
    text(xM0(1)+1.5e7,xM0(2)+1.5e7,xM0(3)-1.e6,txt,'FontSize',14)
    RS = cspice_bodvrd('Sun','RADII',3);
    surf(X*RS(1), Y*RS(2), Z*RS(3),'FaceColor',[0.9290 0.6940 0.1250]), hold off
    txt = {'Sun'};
    text(1.e7,1.e7,1.e6,txt,'FontSize',14)
end

% Results display
tformat = 'YYYY-MON-DD-HR:MN:SC.####::UTC';                                                      
dep_best = cspice_timout(y(1)*TT,tformat); 
arr_best = cspice_timout(y(2)*TT,tformat); 
fprintf(['Best departure date: ',dep_best,'\n']);
fprintf(['Best arrival date: ',arr_best,'\n']);
fprintf('Minimum Δv [km/s]: %f\n',dv);
x0_best = [y(3:5).*LL; y(6:8)];

figure(2)
xE0 = cspice_spkezr(body{1}.name, y(1).*TT, frame, 'NONE', 'SSB');
xM0 = cspice_spkezr(body{2}.name, y(2).*TT, frame, 'NONE', 'SSB');
options_ode = odeset('reltol', 1e-9, 'abstol', 1e-9);
T_Eprop = 688*24*cspice_spd;
[~, xEprop] = ode113(@(t,x) ode_orbit(t,x,mu.*GMGM), [0 T_Eprop], xE0, options_ode); 
T_Mprop = 366*24*cspice_spd;
[~, xMprop] = ode113(@(t,x) ode_orbit(t,x,mu.*GMGM), [0 T_Mprop], xM0, options_ode); 
plot3(xEprop(:,1), xEprop(:,2), xEprop(:,3),'color',[0, 0.5, 0],'LineStyle','-','LineWidth',1), hold on
plot3(xMprop(:,1), xMprop(:,2), xMprop(:,3),'color',[0.8500, 0.3250, 0.0980],'LineStyle','-','LineWidth',1), hold on  
options_ode = odeset('reltol', 1e-8, 'abstol', 1e-8);
[~, y] = ode113(@(t,x) ode_orbit(t,x,data.mu*GMGM), [y(1) y(2)].*TT, x0_best, options_ode);
plot3(y(:,1), y(:,2), y(:,3),'LineStyle','-.','color', 'k','LineWidth',1), hold on
title(sprintf('Best solution found: $\\Delta$v = %f km/s', dv),'Interpreter','latex','FontSize',16)
axis equal
RE = cspice_bodvrd('Earth','RADII',3);
[X, Y, Z] = sphere;
surf(X*RE(1)+xE0(1), Y*RE(2)+xE0(2), Z*RE(3)+xE0(3)), hold on
txt = {'Earth'};
text(xE0(1)-1.5e7,xE0(2)-1.5e7,xE0(3)-1.e6,txt,'FontSize',14)
RM = cspice_bodvrd('Mars','RADII',3);
surf(X*RM(1)+xM0(1), Y*RM(2)+xM0(2), Z*RM(3)+xM0(3)), hold on
txt = {'Mars'};
text(xM0(1)+1.5e7,xM0(2)+1.5e7,xM0(3)-1.e6,txt,'FontSize',14)
RS = cspice_bodvrd('Sun','RADII',3);
surf(X*RS(1), Y*RS(2), Z*RS(3),'FaceColor',[0.9290 0.6940 0.1250]), hold off
txt = {'Sun'};
text(1.e7,1.e7,1.e6,txt,'FontSize',14)
grid on
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')

%% Ex2 
clearvars; close all; clc;

labels = {'Sun';
          'Earth';
          'Moon';};

% Initialize propagation data
bodies = nbody_init(labels);

% select integration frame string (SPICE naming convention)
frame = 'J2000';

TT = 1.e5;
LL = 1.e5;
GMGM = (LL^3)/(TT^2);

data.TT =TT;
data.LL = LL;
data.GMGM = GMGM;
% data.VV = 1; --> not used

mu = bodies{2}.GM./GMGM;

data.labels = labels;
data.frame = frame;
data.mu = mu;
data.h = 200; % altitude of Earth orbit

A = [1,-1,0,zeros(1,12);
     0,1,-1,zeros(1,12)];  
b = [0;
     0];
Aeq = []; % vector velocity along z in Earth frame = 0
beq = [];
tof2 = 20*cspice_spd./TT;
tof1 = 12*cspice_spd./TT;
dep_guess = cspice_str2et('2021-NOV-20-00:00:00')./TT;
med_guess = dep_guess+tof1;
arr_guess = med_guess+tof2;
times = [dep_guess; med_guess; arr_guess];

lb = -Inf(1,15);
ub = Inf(1,15);

xMoon = cspice_spkezr('Moon', times(3).*data.TT, frame, 'NONE', 'Earth');
pos_dep_guess = 6578*(xMoon(1:3)/norm(xMoon(1:3)))./LL;
pos_med_guess = 1500000*(-xMoon(1:3)/norm(xMoon(1:3)))./LL;
vel_dep_guess = 11*(cross([0; 0; 1], pos_dep_guess))/norm(cross([0; 0; 1],pos_dep_guess));
vel_med_guess = .001*(cross([0; 0; 1], pos_med_guess))/norm(cross([0; 0; 1],pos_med_guess));
x0 = [times; pos_dep_guess; pos_med_guess; vel_dep_guess; vel_med_guess];

% Found x0, use this if don't want 11000 function evals which could take 20
% minutes 
% x0 = 1.0e+03 *[6.899272498500432
%                6.918387235856720
%                6.941104491339468
%               -0.000065742845251
%                0.000002250883213
%                                0
%                0.009372134070753
%                0.000896262028757
%               -0.000424619890837
%               -0.000370491756257
%               -0.010963200028622
%                0.000029095143902
%               -0.000337797509790
%                0.000460481454442
%                0.000255666151100];

options = optimoptions('fmincon','Display','iter','Algorithm', ...
    'interior-point','MaxFunctionEvaluations',11000, 'MaxIterations',620);
off = 0; on = 1;
data.turn = off;
data.precision = 1e-10;
[x, dv] = fmincon(@obj_bielliptic, x0, A, b, Aeq, beq, lb, ub, @nonlcons_bielliptic, options, data);
x0 = x;
tof1 = x0(2)-x0(1);
tof2 = x0(3)-x0(2);

figure(2)
xi = cspice_spkezr('Earth', x(1)*data.TT, frame, 'NONE', 'Earth');
xMoon_i = cspice_spkezr('Moon', x(1)*data.TT, frame, 'NONE', 'Earth');
xMoon_m = cspice_spkezr('Moon', x(2)*data.TT, frame, 'NONE', 'Earth');
xMoon_f = cspice_spkezr('Moon', x(3)*data.TT, frame, 'NONE', 'Earth');

options_ode = odeset('reltol', 1e-12, 'abstol', 1e-12);
x1 = x(4:6);       x1(4:6) = x(10:12);
xm = x(7:9);       xm(4:6) = x(13:15);
[t1, x_tran1] = ode113(@(t,x) nbody_rhs(t,x,bodies,frame,data), [x(1) x(2)], x1, options_ode);
[t2, x_tran2] = ode113(@(t,x) nbody_rhs(t,x,bodies,frame,data), [x(2) x(3)], xm, options_ode); 
t = [t1; t2];
x_Moon = zeros((length(x_tran1)+length(x_tran2)),6);
for i = 1:(length(x_tran1)+length(x_tran2))
    x_Moon(i,:) = cspice_spkezr('Moon', t(i)*data.TT, frame, 'NONE', 'Earth');
end
tformat = 'YYYY-MON-DD-HR:MN:SC.####::UTC';
dep_best = cspice_timout(x(1)*data.TT,tformat); 
arr_best = cspice_timout(x(3)*data.TT,tformat); 

figure(2)
plot3(x_Moon(:,1), x_Moon(:,2), x_Moon(:,3),'Color','k','LineStyle','-.','LineWidth',1.2)
xL2 = (xMoon_f(1:3)+60000*xMoon_f(1:3)/norm(xMoon_f(1:3),2)); 
scatter3(xL2(1), xL2(2),  xL2(3), 7.e1,'Marker','hexagram','MarkerFaceColor','flat','DisplayName','L2','color',[0.4940, 0.1840, 0.5560],'Displayname', 'L2'), hold on
scatter3(xm(1)*data.LL, xm(2)*data.LL,  xm(3)*data.LL, 7.e1,'Marker','*','DisplayName','Manoeuvre','color',[0.9290, 0.6940, 0.1250],'Displayname', 'Manoeuvre point','LineWidth',1.2), hold on
scatter3(x_Moon(1,1), x_Moon(1,2), x_Moon(1,3),5.e1,'LineWidth',1.3, 'Marker','o','color', [0, 0, 1],'Displayname',['Moon at ', dep_best]), hold on
scatter3(x_Moon(end,1), x_Moon(end,2), x_Moon(end,3), 7.e1,'LineWidth',1.3, 'Marker', 'x','color', [1, 0, 0],'Displayname', ['Moon at ', arr_best]), hold on
plot3(x_Moon(:,1), x_Moon(:,2), x_Moon(:,3),'Color',[0.5 0.5 0.5],'LineStyle','-.','LineWidth',1.2)
title(sprintf('Best solution found: $\\Delta$v = %f km/s', dv),'Interpreter','latex')
axis equal
grid on 
legend show, legend('AutoUpdate','off')
plot3(x_tran1(:,1)*data.LL, x_tran1(:,2)*data.LL, x_tran1(:,3)*data.LL,'k','LineWidth',1.2), hold on
plot3(x_tran2(:,1)*data.LL, x_tran2(:,2)*data.LL, x_tran2(:,3)*data.LL,'k','LineWidth',1.2), hold on 
RE = cspice_bodvrd('Earth','RADII',3);
[X, Y, Z] = sphere;
surf(X*RE(1), Y*RE(2), Z*RE(3)), hold on
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')

% Results display                                                      
fprintf(['Best departure date: ',dep_best,'\n']);
fprintf(['Best arrival date: ',arr_best,'\n']);
TOF1 = (x0(2)-x0(1))*data.TT/cspice_spd;
TOF2 = (x0(3)-x0(2))*data.TT/cspice_spd;
fprintf(['TOF1: ',num2str(TOF1),'\n']);
fprintf(['TOF2: ',num2str(TOF2),'\n']);
fprintf('Minimum Δv [km/s]: %f\n',dv);

N = 20;
prompt = '\n Part 3 may take an hour. Continue? Y/N [N]: ';
str = input(prompt,'s');
if isempty(str) || str == 'N' || str == 'n'
    str = 'N';
end
if str== 'Y' || str== 'y'
   options = optimoptions('fmincon','Display','iter','Algorithm', ...
    'interior-point','MaxFunctionEvaluations',1200, 'MaxIterations',80, ...
    'StepTolerance',1.e-4,'OptimalityTolerance',1.e-4,'FunctionTolerance', ...
    1.e-4,'ConstraintTolerance',1.e-4);
    Dv = zeros(N,1);
    exitflag = zeros(N,1);
    data.precision = 1.e-8;
    DEP = linspace(x0(1), x0(1)+365*cspice_spd/data.TT, N);
    for i = 1:N
        med_guess = DEP(i)+tof1;
        arr_guess = med_guess+tof2;
        times = [DEP(i); med_guess; arr_guess];
        xMoon = cspice_spkezr('Moon', times(3).*data.TT, frame, 'NONE', 'Earth');
        pos_dep_guess = 6578*(xMoon(1:3)/norm(xMoon(1:3)))./LL;
        pos_med_guess = 1500000*(-xMoon(1:3)/norm(xMoon(1:3)))./LL;
        vel_dep_guess = 11*(cross([0; 0; 1], pos_dep_guess))/norm(cross([0; 0; 1],pos_dep_guess));
        vel_med_guess = .001*(cross([0; 0; 1], pos_med_guess))/norm(cross([0; 0; 1],pos_med_guess));
        x0 = [times; pos_dep_guess; pos_med_guess; vel_dep_guess; vel_med_guess];
        dep_guess_max = DEP(i)+6000/data.TT;
        dep_guess_min = DEP(i)-6000/data.TT;
        lb = [dep_guess_min, -Inf(1,14)];
        ub = [dep_guess_max, Inf(1,14)];
        [x, dv, flag] = fmincon(@obj_bielliptic, x0, A, b, Aeq, beq, lb, ub, @nonlcons_bielliptic, options, data);
        Dv(i) = dv;
        DEP(i) = x(1);
        exitflag(i) = flag;
    end
    DEP = linspace(x0(1), x0(1)+365*cspice_spd/data.TT, N);
    dep_UTC = cspice_timout(DEP*data.TT,tformat);
    dn = datenum(dep_UTC, 'yyyy-mmm-dd');
    figure()
    plot(dn, Dv,'LineWidth',1.1,'DisplayName','Output'), hold on
    datetick('x', 'yyyy-mmm-dd','keepticks')
    axis tight
    title('Earth-L2 (EML2) transfer','FontSize', 14)
    grid on
    xlabel('Departure date','FontSize', 14)
    ax=gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    ylabel('Δv [km/s]','FontSize', 14)
    hold on
    plot(dn, movmean(Dv, [0 4]), 'LineStyle','--','color', [0.8500, 0.3250, 0.0980],'LineWidth',1.1,'DisplayName','Movemean with 4 points')
    plot(dn, movmean(Dv, [0 22]), 'LineStyle','-.','color', 'black','LineWidth',1.1,'DisplayName','Movemean with 22 points')
    legend show

end


%% Ex 3
clearvars; close all; clc;

% Initialize propagation data
bodies = nbody_init({'Earth'});

frame = 'J2000';

mu = bodies{1}.GM;

% Data
r0 = [0; -29597.43; 0];    v0 = [1.8349; 0.0002; 3.1783];    m0 = 735;   
rf = [0; -29617.43; 0];    vf = [1.8371; 0.0002; 3.1755];
Tmax = 100.e-6;        % kN
Isp = 3000;            % s
g0 = 9.8*1.e-3;        % km/s^2

lambda0 = guess_lambda;
data.Tmax = Tmax;
data.Isp = Isp;
data.g0 = g0;
data.xf = [rf; vf];
data.x0 = [r0; v0; m0];
data.mu = mu;

% Non-dimensionalization
data.TT = 1.e0;
data.LL = 1.e0;
data.MM = 1.e0;
data.VV = data.LL/data.TT;
data.acc = data.LL/data.TT^2;
data.GMGM = data.LL^3/data.TT^2;

data.Tmax = data.Tmax/(data.MM*data.acc);
data.Isp = data.Isp/data.TT;
data.g0 = data.g0/data.acc;
data.mu = data.mu/data.GMGM;
data.x0(1:3) = data.x0(1:3)/data.LL;
data.x0(4:6) = data.x0(4:6)/data.VV;
data.x0(7) = data.x0(7)/data.MM;
data.xf(1:3) = data.xf(1:3)/data.LL;
data.xf(4:6) = data.xf(4:6)/data.VV;

T = 2*pi*sqrt(norm(rf,2)^3/mu); % guess since the semi-major axis is reduced

data.lR = 1/data.LL;
data.lV = 1/data.VV;
data.lM = 1/data.MM;

% Solver
y0 = [T/data.TT; lambda0(1:3)/data.lR; lambda0(4:6)/data.lV; lambda0(7)/data.lM];

% y0 to have the same results as in the report:
% y0 = 1.0e+04 *[ 5.072616616422452
%                0.000444879184734
%               -0.000700269115044
%                0.000319210505817
%                0.518594937696488
%                0.972974554493608
%                0.648991489202271
%                0.009001652876762];

options = optimoptions(@fsolve,'Algorithm', 'Levenberg-Marquardt', 'Display', 'iter','Maxiter', 1000,'FunctionTolerance', 1e-6, 'StepTolerance',1e-6);

zero = fsolve(@low_thrust_zero, y0, options, data); 

% Display results
figure(2)
options_ode = odeset('reltol', 1.e-13, 'abstol', 1e-16);
[~, y] = ode113(@(t,x) lt_integration(t, x, data), [0 zero(1)], [data.x0; zero(2:8)], options_ode);
plot3(y(:,1)*data.LL, y(:,2)*data.LL, y(:,3)*data.LL,'k'), hold on
title(sprintf('Minimum time = %f s', zero(1)*data.TT),'Interpreter','latex')
scatter3(y(1,1)*data.LL, y(1,2)*data.LL, y(1,3)*data.LL,3.e1,'b','Marker','o','LineWidth',1.3), hold on  
scatter3(y(end,1)*data.LL, y(end,2)*data.LL, y(end,3)*data.LL, 8.e1, 'r','Marker','x','LineWidth',1.3), hold on
axis equal
txt = {'Start'};
text(y(1,1)*data.LL-4.e1,y(1,2)*data.LL+2e1,y(1,3)*data.LL+1e1,txt,'FontSize',13)
txt = {'Target'};
text(y(1,1)*data.LL+2.e1,y(1,2)*data.LL-1e1,y(1,3)*data.LL-1e1,txt,'FontSize',13)
xlim([-50 50])
ylim(1.e4*[-2.969 -2.95])
zlim([-150 100])
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
grid on

zero(1) = zero(1)*data.TT;
zero(2:4) = zero(2:4)*data.lR;
zero(5:7) = zero(5:7)*data.lV;
zero(8) = zero(8)*data.lM;
fprintf(['The solution of the following problem is:\n \n λ = \n [', ...
    num2str(zero(2)), '\n', num2str(zero(3)), '\n', ...
    num2str(zero(4)), '\n', num2str(zero(5)), '\n', num2str(zero(6)), '\n', ...
    num2str(zero(7)), '\n', num2str(zero(8)),']\n']); 
fprintf('\n and \n \n tf = %f\n \n', zero(1));

%% functions 
% Ex 2
function [dv] = obj_bielliptic(x, data)
%
%     obj_bielliptic - Δv objective for a bielliptic trajectory 
%
%     DESCRIPTION:
%       Function that gives as output the objective for a Δv minimization
%       in a bielliptic trajectory case from Earth orbit to Lagrangian 
%       point L2
%
%     PROTOTYPE:
%       [dv] = obj_bielliptic(x, data)
%          
%     INPUT:
%       x[15,1]        Minimizing variables
%       data           Struct with fields:
%                      - data.labels: e.g. {'Earth'}
%                      - data.TT: time non-dimensionalization
%                      - data.frame: 'J2000' 
%                      - data.LL: space non-dimensionalization
%                      - data.GMGM: gravitational constant
%                                   non-dimensionalization
%                      - data.h: altitude of Earth orbit
%                      - data.turn: display every iteration (advice just to
%                                   comprehend trajctory since makes the
%                                   program slower), = 0(off)/1(on)
%                      - data.precision: reltol and abstol required
%     
%     OUTPUT:
%       dv              Δv required for the bielliptic selected%                      
%     
%     CALLED FUNCTIONS:
%      nbody_init, cspice_bodvrd, cspice_spkezr, nbody_rhs
% 
%     LAST UPDATED:
%      8/11/2021
%
%     CREATED BY:
%      Bardazzi Niccolò

ti = x(1);         tm = x(2);         tf = x(3);
x1 = x(4:6);       x1(4:6) = x(10:12);
xm = x(7:9);       xm(4:6) = x(13:15);
bodies = nbody_init(data.labels);
precision = data.precision;

% L2
xMoon = cspice_spkezr('Moon', tf.*data.TT, data.frame, 'NONE', 'Earth');
posL2 = (xMoon(1:3)+60000*xMoon(1:3)/norm(xMoon(1:3),2))./data.LL; 
velL2 = xMoon(4:6)*norm(posL2)/norm(xMoon(1:3))*data.LL;
xL2 = [posL2; velL2];

% Initial position in orbit
RE = cspice_bodvrd('Earth','RADII',3);
v_orb = sqrt(data.mu.*data.GMGM/(RE(1)+data.h));

% Inertial integration
options = odeset('reltol', precision, 'abstol', precision);
[~, x_transfer1] = ode113(@(t,x) nbody_rhs(t,x,bodies,'J2000',data), [ti tm], x1, options);                                            
dv1 = x1(4:6)-v_orb*cross([0; 0; 1], x1(1:3))/norm(cross([0; 0; 1], x1(1:3)),2);
[~, x_transfer2] = ode113(@(t,x) nbody_rhs(t,x,bodies,'J2000',data), [tm tf], xm, options);
dv2 = xm(4:6)-x_transfer1(end,4:6)';
dv3 = xL2(4:6)-x_transfer2(end,4:6)';
dv = norm(dv1,2)+norm(dv2,2)+norm(dv3,2);

if data.turn == 1
    x_transfer1(:,1:3) = x_transfer1(:,1:3).*data.LL;
    x_transfer2(:,1:3) = x_transfer2(:,1:3).*data.LL;
    figure(1)
    scatter3(0,0,0,'Marker','diamond', 'DisplayName', 'Earth'), hold on
    plot3(x_transfer1(:,1), x_transfer1(:,2), x_transfer1(:,3), 'DisplayName', 'Transfer 1','LineStyle','-.'), hold on
    plot3(x_transfer2(:,1), x_transfer2(:,2), x_transfer2(:,3), 'DisplayName', 'Transfer 2'), hold on
    scatter3(xm(1).*data.LL, xm(2).*data.LL, xm(3).*data.LL, 'Marker', '*', 'DisplayName', 'Manoeuvre point')
    scatter3(xL2(1).*data.LL, xL2(2).*data.LL,  xL2(3).*data.LL,'Marker','hexagram','DisplayName', 'L2'), hold on
    scatter3(xMoon(1), xMoon(2), xMoon(3),'k', 'DisplayName', 'Moon'), hold on 
    xSun = cspice_spkezr('Sun', tf.*data.TT, data.frame, 'NONE', 'Earth');
    quiver3(0,0,0,xSun(1)/20, xSun(2)/20, xSun(3)/20,'DisplayName', 'Sun direction'), hold on
    legend show, legend('AutoUpdate','off')
    axis equal
end

end

function [c, ceq] = nonlcons_bielliptic(x, data)
%
%     nonlcons_bielliptic - Non-linear constraints for a bielliptic trajectory 
%
%     DESCRIPTION:
%       Function that gives as output the non-linear constraint, equality
%       and disequality, for a bielliptic trajectory from Earth orbit to
%       Lagrangian point L2
%
%     PROTOTYPE:
%       [c, ceq] = nonlcons_bielliptic(x, data)
%          
%     INPUT:
%       x[15,1]        Minimizing variables
%       data           Struct with fields:
%                      - data.labels: e.g. {'Earth'}
%                      - data.TT: time non-dimensionalization
%                      - data.frame: 'J2000' 
%                      - data.LL: space non-dimensionalization
%                      - data.GMGM: gravitational constant
%                                   non-dimensionalization
%                      - data.h: altitude of Earth orbit
%                      - data.precision: reltol and abstol required
%     
%     OUTPUT:
%       c              Minimum altitude 500000 km
%       ceq            [Match between branches at manoeuvering point; 
%                       match at the arrival point with L2;
%                       match at the parking orbit in altitude norm;
%                       match at the parking orbit as equatorial orbit;]
%                      
%     
%     CALLED FUNCTIONS:
%      nbody_init, cspice_bodvrd, cspice_spkezr, nbody_rhs
% 
%     LAST UPDATED:
%      8/11/2021
%
%     CREATED BY:
%      Bardazzi Niccolò

    ti = x(1);     tm = x(2);     tf = x(3);
    x1 = x(4:6);       x1(4:6) = x(10:12);
    xm = x(7:9);       xm(4:6) = x(13:15); 
    bodies = nbody_init(data.labels);
    precision = data.precision;
    
    % L2
    xMoon = cspice_spkezr('Moon', tf.*data.TT, data.frame, 'NONE', 'Earth');
    posL2 = (xMoon(1:3)+60000*xMoon(1:3)/norm(xMoon(1:3),2))./data.LL; 
    velL2 = xMoon(4:6)*norm(posL2)/norm(xMoon(1:3))*data.LL;
    xL2 = [posL2; velL2];

    % norm of position and z component of position
    RE = cspice_bodvrd('Earth','RADII',3);

    % Inertial integration
    options = odeset('reltol', precision, 'abstol', precision);
    [~, x_transfer1] = ode113(@(t,x) nbody_rhs(t,x,bodies,'J2000',data), [ti tm], x1, options);
    [~, x_transfer2] = ode113(@(t,x) nbody_rhs(t,x,bodies,'J2000',data), [tm tf], xm, options);
    ceq = [xm(1:3)-x_transfer1(end,1:3)'; xL2(1:3)-x_transfer2(end,1:3)'; norm(x1(1:3),2)-(data.h+RE(1))./data.LL; x1(3)];
    c = 500000./data.LL-norm(x_transfer1(end,1:3),2);
end

function [dx] = nbody_rhs(t, x, bodies, frame, data)

%NBODY_RHS Evaluates the right-hand-side of a N-body propagator
% Evaluates the right-hand-side of a newtonian N-body propagator.
% The integration centre is the Solar-System-Barycentre (SSB) and only
% Newtonian gravitational accelerations are considered.
%
%
% Author
% Name: ALESSANDRO
% Surname: MORSELLI
% Research group: DART
% Department: DAER
% University: Politecnico di Milano
% Creation: 26/09/2021
% Contact: alessandro.morselli@polimi.it
% Copyright: (c) 2021 A. Morselli, Politecnico di Milano.
% All rights reserved.
%
%
% Notes:
% This material was prepared to support the course 'Satellite Guidance
% and Navigation', AY 2021/2022.
%
%
% Inputs:
% t : [1,1] ephemeris time (ET SPICE), seconds past J2000 (TDB)
% x : [6,1] cartesian state vector wrt Solar-System-Barycentre
% bodies : [1,6] cell-array created with function nbody_init
%
% Outputs:
% dxdt : [6,1] RHS, newtonian gravitational acceleration only
%
% Prerequisites:
% - MICE (Matlab SPICE)
% - Populated kernel pool (SPK, LSK, PCK kernels)
%
%     LAST UPDATED:
%      8/11/2021
%
%     BY:
%      Bardazzi Niccolò
%
%     NEW PROTOTYPE:
%       [dx] = nbody_rhs(t, x, bodies, frame, data)
%          
%     ADDITIONAL INPUT:
%       data           Struct with fields:
%                      - data.TT: time non-dimensionalization
%                      - data.LL: space non-dimensionalization
%                      - data.GMGM: gravitational constant
%                                   non-dimensionalization
% 

if not( strcmpi(frame, 'ECLIPJ2000') || strcmpi(frame, 'J2000') )
msg = 'Invalid integration reference frame, select either J2000 or ECLIPJ2000';
error(msg);
end

% Initialize right-hand-side
dx = zeros(6,1);
% Position derivative is object's velocity
dx(1:3) = x(4:6);

% Extract the object position from state x
rr_ssb_obj = x(1:3); 
for i=1:length(bodies)

% Retrieve position and velocity of i-th celestial body wrt Solar
% System Barycentre in inertial frame
rv_ssb_body = cspice_spkezr(bodies{i}.name, t*data.TT, frame, 'NONE', 'Earth');

% Extract object position wrt. i-th celestial body
rr_body_obj = rr_ssb_obj.*data.LL - rv_ssb_body(1:3);

% Compute square distance and distance
dist2 = dot(rr_body_obj, rr_body_obj);
dist = sqrt(dist2);

xS = cspice_spkezr('SUN', t*data.TT, frame, 'NONE','Earth');
xM = cspice_spkezr('MOON', t*data.TT, frame, 'NONE','Earth');
rS = xS(1:3); rS3 = dot(rS, rS)*sqrt(dot(rS,rS));
rM = xM(1:3); rM3 = dot(rM, rM)*sqrt(dot(rM,rM));
rE = [0;0;0]; rE3 = 1;
r = [rS rE rM];
r3 = [rS3; rE3; rM3];

% Compute the gravitational acceleration using Newton's law
aa_grav = - bodies{i}.GM * (rr_body_obj /(dist*dist2) + r(:,i)./r3(i));

% Sum up acceleration to right-hand-side
dx(4:6) = dx(4:6) + aa_grav.*data.GMGM;

end

end

% Ex 3
function lambda0 = guess_lambda 
%
%     guess_lambda - Initial guess of lambda for OCP 
%
%     DESCRIPTION:
%       Function that gives an initial estimation of the Lagrangian 
%       multipliers for on Earth OCP problem
%
%     PROTOTYPE:
%       lambda0 = guess_lambda 
%          
%     INPUT:
%       -
%     
%     OUTPUT:
%       lambda0[7,1]   Vector full of random numbers such that:
%                      - -10<lambda(1:3,1)<10
%                      - -10^4<lambda(4:6,1)<10^4
%                      - 50<lambda(7,1)<100
%     
%     CALLED FUNCTIONS:
%       -
% 
%     LAST UPDATED:
%      8/11/2021
%
%     CREATED BY:
%      Bardazzi Niccolò

pos_min = -10;
pos_max = 10;
lambda_pos = pos_min+rand(3,1)*(pos_max-pos_min);

vel_min = -1.e4;
vel_max = 1.e4;
lambda_vel = vel_min+rand(3,1)*(vel_max-vel_min);

mass_min = 50;
mass_max = 100;
lambda_mass = mass_min+rand(1)*(mass_max-mass_min);

lambda0 = [lambda_pos; lambda_vel; lambda_mass];
end

function fx = low_thrust_zero(z, data)
%
%     low_thrust_zero - Low thrust function for moving from a position to 
%                       an other
%
%     DESCRIPTION:
%       Function to be given to a vectorial zero finding solver (fsolve in
%       MATLAB) to recover, with a shooting method, the values of the 
%       Lagrangian multipliers 
%
%     PROTOTYPE:
%       fx = low_thrust_zero(z, data)
%          
%     INPUT:
%       z[8,1]         Vector of zeros:
%                      - time required, z(1,1)
%                      - Lagrangian multipliers (2:8,1)
%       data           Struct with fields:
%                      - data.xf: final state
%                      - data.x0: initial state
%                      - data.Tmax: maximum thrust
%                      - data.Isp: specific impulse
%                      - data.mu: gravitational constant 
%                      - data.g0: 9.8 m/s^2
%     
%     OUTPUT:
%       fx             Vector to be minimized: [x(t_f)-xf;dλmdt(t_f); H]
%     
%     CALLED FUNCTIONS:
%       -
% 
%     LAST UPDATED:
%      8/11/2021
%
%     CREATED BY:
%      Bardazzi Niccolò

    xf = data.xf;
    x0 = data.x0;
    
    options = odeset('reltol', 1e-12, 'abstol', 1e-12);
    [~, y] = ode113(@(t,x) lt_integration(t, x, data), [0 z(1)], [x0; z(2:8)], options);
    fx = [(y(end,1:3)'-xf(1:3)); y(end,4:6)'-xf(4:6); y(end,14)];

end

function dxdt = lt_integration(~, x, data)
%
%     lt_integration - Low thrust for time optimal problem system to be 
%                      integrated
%
%     DESCRIPTION:
%       Function to be provided to the ODE integrator for the integration
%       of an orbit with low thrust propulsion for a time optimal problem
%
%     PROTOTYPE:
%       dxdt = lt_integration(~, x, data)
%          
%     INPUT:
%       t              Integration time 
%       y              Integration vector (the one with initial conditions)
%       data           Struct with fields:
%                      - data.Tmax: maximum thrust
%                      - data.Isp: specific impulse
%                      - data.mu: gravitational constant 
%                      - data.g0: 9.8 m/s^2
%     
%     OUTPUT:
%       dxdt           Vector to be integrated
%     
%     CALLED FUNCTIONS:
%       -
% 
%     LAST UPDATED:
%      8/11/2021
%
%     CREATED BY:
%      Bardazzi Niccolò

u = 1;
Tmax = data.Tmax;
Isp = data.Isp;
g0 = data.g0;
mu = data.mu;
dxdt = zeros(14,1);

% r = x(1:3);            v = x(4:6);              m = x(7);
% lambdar = x(8:10);     lambdav = x(11:13);      % lambdam = x(14);      
% dxdt(1:3) = v;
% dxdt(4:6) = -mu/norm(r,2)^3*r-u*Tmax/m*lambdav/norm(lambdav,2);
% dxdt(7) = -u*Tmax/(Isp*g0);
% dxdt(8:10) = -3*mu/norm(r,2)^5*r'*lambdav*r+mu/norm(r,2)^3*lambdav;
% dxdt(11:13) = -lambdar;
% dxdt(14) = -u*norm(lambdav,2)*Tmax/m^2;

r = x(1:3); v = x(4:6); m = x(7);
lambda_r = x(8:10); lambda_v = x(11:13);
dxdt(1:3) = v;
dxdt(4:6) = -mu/norm(r)^3.*r-u*Tmax/m.*lambda_v./norm(lambda_v);
dxdt(7) = -u*Tmax/(Isp*g0);
dxdt(8:10) = -3*mu/norm(r)^5*dot(r, lambda_v).*r+mu/norm(r)^3.*lambda_v;
dxdt(11:13) = -lambda_r;
dxdt(14) = -u*norm(lambda_v)*Tmax/m^2;

end

% Ex 1
function [dv] = obj(x, data)
%
%     obj - objective function aimed at minimizing Δv
%
%     DESCRIPTION:
%       Objective function to be provided to an optimization solver in
%       order to minimize the Δv of a trajectory with patched conic method
%
%     PROTOTYPE:
%        [dv] = obj(x, data)
%          
%     INPUT:
%       x[8,1]         Variables of the minimization probelem: 
%                      - x(1): scaled depature time
%                      - x(2): scaled arrival time
%                      - x(3:5): scaled departure position
%                      - x(6:8): scaled departure velocity
%       data           Struct with fields:
%                      - data.TT: scaling in time
%                      - body{}: struct with the 2 planets name, example:
%                        {'Earth'}, {'Mars'}
%                      - data.LL: scaling in position
%                      - data.frame: 'J2000' 
%                      - data.mu: gravitational constant of the main
%                        attractor
%                      - data.turn: display in plotting trajectories
%     
%     OUTPUT:
%       dv             minimum Δv of the trajectory
%     
%     CALLED FUNCTIONS:
%      ode_orbit, cspice_spkezr
% 
%     LAST UPDATED:
%      7/11/2021
%
%     CREATED BY:
%      Bardazzi Niccolò

x1 = x(3:8);
ti = x(1)*data.TT;                tf = x(2)*data.TT;   
body{1} = data.body{1};           body{2} = data.body{2};
frame = data.frame;
mu = data.mu;

xi = cspice_spkezr(body{1}.name, ti, frame, 'NONE', 'SSB');
xf = cspice_spkezr(body{2}.name, tf, frame, 'NONE', 'SSB');

xi(1:3) = xi(1:3)./data.LL;
xf(1:3) = xf(1:3)./data.LL;

options_ode = odeset('reltol', 1e-9, 'abstol', 1e-9);
[~, y] = ode113(@(t,x) ode_orbit(t, x, mu), [ti tf]./data.TT, x1, options_ode);
x1 = y(1,1:6)';            x2 = y(end,1:6)';
dv1 = x1(4:6)-xi(4:6);
dv2 = xf(4:6)-x2(4:6);
% objective 
dv = norm(dv1,2) + norm(dv2,2);

y(:,1:3) = y(:,1:3).*data.LL;
if data.turn == 1
    figure(1)
    scatter3(y(end,1), y(end,2), y(end,3),'Marker','hexagram','MarkerFaceColor','flat'), hold on
    scatter3(y(1,1), y(1,2), y(1,3),'Marker','o'), hold on
    plot3(y(:,1), y(:,2), y(:,3)), hold on
    axis equal
    grid on
end
end

function [c, ceq] = nonlcons(x, data)
%
%     nonlcons - nonlinear constraint of the problem aimed at minimizing Δv
%
%     DESCRIPTION:
%       Objective function to be provided to an optimization solver in
%       order to minimize the Δv of a trajectory with patched conic method
%
%     PROTOTYPE:
%        [dv] = obj(x, data)
%          
%     INPUT:
%       x[8,1]         Variables of the minimization probelem: 
%                      - x(1): scaled depature time
%                      - x(2): scaled arrival time
%                      - x(3:5): scaled departure position
%                      - x(6:8): scaled departure velocity
%       data           Struct with fields:
%                      - data.TT: scaling in time
%                      - body{}: struct with the 2 planets name, example:
%                        {'Earth'}, {'Mars'}
%                      - data.LL: scaling in position
%                      - data.frame: 'J2000' 
%                      - data.mu: gravitational constant of the main
%                        attractor
%                      - data.turn: option in plotting trajectories
%     
%     OUTPUT:
%       c[]            Linear constraints, empty
%       ceq[2,1]       Non-linear constraints aimed at matching initial
%                      position and final position
%     
%     CALLED FUNCTIONS:
%      ode_orbit, cspice_spkezr
% 
%     LAST UPDATED:
%      7/11/2021
%
%     CREATED BY:
%      Bardazzi Niccolò

x1 = x(3:8);
ti = x(1);                tf = x(2);   
body{1} = data.body{1};   body{2} = data.body{2};
frame = data.frame;
mu = data.mu;

xi = cspice_spkezr(body{1}.name, ti.*data.TT, frame, 'NONE', 'SSB');
xf = cspice_spkezr(body{2}.name, tf.*data.TT, frame, 'NONE', 'SSB');
xi(1:3) = xi(1:3)./data.LL;
xf(1:3) = xf(1:3)./data.LL;

options_ode = odeset('reltol', 1e-9, 'abstol', 1e-9);
[~, y] = ode113(@(t,x) ode_orbit(t,x, mu), [ti tf], x1, options_ode);
x1 = y(1,1:6)';            x2 = y(end,1:6)';

ceq = [x1(1:3)-xi(1:3); x2(1:3)-xf(1:3)];
c = [];
end

function [x1, x2, count] = NM_STM(x1_guess, x_target, mu, tspan, tol, options_ode)
%
%     NM_STM - zero finding algorithm for Lambert's problem solution
%
%     DESCRIPTION:
%       Newton Method algorithm based on variational approach, thus on STM 
%       (State Transition Matrix) aimed in the zero finding of the final 
%       position
%
%     PROTOTYPE:
%        [x1, x2, count] = NM_STM(x1_guess, x_target, mu, tspan, tol, options_ode)
%          
%     INPUT:
%       x1_guess[6,1]  Initial state: [correct position; guess in velocity] 
%       x_target[3,1]  Target position
%       mu             Gravitational constant
%       tspan[1,2]     Time span 
%       tol            Tolerance to be satisfied
%       options_ode    ODE options (MATLAB ode)
%     
%     OUTPUT:
%       x1             Correct initial state
%       x2             Final state
%       count          N° of iterations
%     
%     CALLED FUNCTIONS:
%      ode_STM
% 
%     LAST UPDATED:
%      6/11/2021
%
%     CREATED BY:
%      Bardazzi Niccolò

drf = ones(3,1);
count = 0;
maxcount = 100;
while norm(drf,2)>tol && count<maxcount
    I = eye(6);
    yPhi0 = [x1_guess; I(:)];
    [~, y] = ode113(@(t,x) ode_STM(t,x, mu), tspan, yPhi0, options_ode);
    x2 = y(end,1:6);
    STM = reshape(y(end, 7:42), [6, 6]);
    % dr0 = 0! --> only Phi_rv is needed
    STMrv = STM(1:3,4:6); 
    drf = x2(1:3)'-x_target;
    dv0 = STMrv\drf;
    x1_guess = [y(1,1:3)'; y(1,4:6)'-dv0];
    count = count+1;
end
x2 = x2';
x1 = x1_guess;
if count == maxcount % no convergence
    x1 = NaN(6,1);
    disp('No solution found')
end
end

function dy = ode_STM(~, yPhi, mu)
%
%     ode_STM - Function of integration
%
%     DESCRIPTION:
%       Function to be provided to the ODE integrator for STM (State
%       Transition Matrix) problem
%
%     PROTOTYPE:
%        dy = ode_STM(~, yPhi, mu)
%          
%     INPUT:
%       t              Initial state: [correct position; guess in velocity] 
%       yPhi           Integration vector (the one with initial conditions)
%       mu             Gravitational constant
%     
%     OUTPUT:
%       dy             Vector to be integrated
%     
%     CALLED FUNCTIONS:
%       -
% 
%     LAST UPDATED:
%      6/11/2021
%
%     CREATED BY:
%      Bardazzi Niccolò

    r = yPhi(1:3);
    v = yPhi(4:6);
    Phi = yPhi(7:42);
    Phi = reshape(Phi, [6,6]);
    A = A_mat(r, mu);
    Phi_dot = A*Phi;
    x_dot = zeros(6,1);
    x_dot(1:3,1) = v;
    x_dot(4:6,1) = -mu*r/norm(r,2).^3;
    dy = [x_dot; Phi_dot(:)];
end

function [x1, x2, count] = NM_FD(x1_guess, x_target, mu, tspan, tol, options_ode)
%
%     NM_FD - zero finding algorithm for Lambert's problem solution
%
%     DESCRIPTION:
%       Newton Method algorithm based on forward Finite Differences aimed 
%       in the zero finding of the final position
%
%     PROTOTYPE:
%        [x1, x2, count] = NM_FD(x1_guess, x_target, mu, tspan, tol, options_ode)
%          
%     INPUT:
%       x1_guess[6,1]  Initial state: [correct position; guess in velocity] 
%       x_target[3,1]  Target position
%       mu             Gravitational constant
%       tspan[1,2]     Time span 
%       tol            Tolerance to be satisfied
%       options_ode    ODE options (MATLAB ode)
%     
%     OUTPUT:
%       x1             Correct initial state
%       x2             Final state
%       count          N° of iterations (max 100)
%     
%     CALLED FUNCTIONS:
%      ode_orbit
% 
%     LAST UPDATED:
%      6/11/2021
%
%     CREATED BY:
%      Bardazzi Niccolò

drf = ones(3,1);
count = 0;
maxcount = 100;
while norm(drf,2)>tol && count<maxcount
    [~, y] = ode113(@(t,x) ode_orbit(t,x,mu), tspan, x1_guess, options_ode);
    x2 = y(end,1:6)';
    STM = diffjac(x1_guess, x2, tspan, mu, options_ode);
    STMrv = STM(1:3,4:6);
    drf = x2(1:3)-x_target;
    dv0 = STMrv\drf;
    x1_guess = [y(1,1:3)'; y(1,4:6)'-dv0]; 
    count = count+1;
end
x1 = x1_guess;
if count>maxcount % convergence not reached
    x1 = NaN;
    disp('No solution found')
end
end

function A = A_mat(r, mu)
%
%     A_mat - Matrix of the exact dyamics of the 2-body problem
%
%     DESCRIPTION:
%       Matrix generation useful for variational approach which provides
%       the exact dynamics of the 2-body problem wrt state space model
%
%     PROTOTYPE:
%        A = A_mat(r, mu)
%          
%     INPUT:
%       r[3,1]         Position vector
%       mu             Gravitational constant
%     
%     OUTPUT:
%       A[6,6]         Matrix of dynamics
%     
%     CALLED FUNCTIONS:
%       -
% 
%     LAST UPDATED:
%      6/11/2021
%
%     CREATED BY:
%      Bardazzi Niccolò

a = 3*mu/norm(r, 2)^5*[r(1)^2, r(1)*r(2), r(1)*r(3); r(1)*r(2), r(2)^2, r(2)*r(3); r(1)*r(3), r(2)*r(3), r(3)^2]-mu/norm(r, 2)^3*eye(3);
A = [ zeros(3)   eye(3)
         a      zeros(3) ];
end

function Jf_FD = diffjac(x0, Phi0, tspan, mu, options_ode)
%
%     diffjac - jacobian of the flux computed with forward finite
%               differences
%
%     DESCRIPTION:
%       Jacobian of the flux of a 2-body problem computed with forward 
%       finite differences useful for Newton's method 
%
%     PROTOTYPE:
%        Jf_FD = diffjac(x0, Phi0, tspan, mu, options_ode)
%          
%     INPUT:
%       x0[6,1]        Initial state: [correct position; guess in velocity] 
%       Phi0[6,1]      Final state
%       tspan[1,2]     Time span 
%       mu             Gravitational constant
%       options_ode    ODE options (MATLAB ode)
%     
%     OUTPUT:
%       x1             Correct initial state
%       x2             Final state
%       count          N° of iterations (max 100)
%     
%     CALLED FUNCTIONS:
%      ode_orbit
% 
%     LAST UPDATED:
%      6/11/2021
%
%     CREATED BY:
%      Bardazzi Niccolò

theta = max(sqrt(eps), sqrt(eps)*abs(x0));
Jf_FD = zeros(length(x0), length(x0));
for j = 1:length(x0)
    pert = zeros(6,1);
    pert(j) = theta(j);
    [~, xx_pert] = ode113(@(t,x) ode_orbit(t,x,mu), tspan, x0+pert, options_ode);
    Phi_pert = xx_pert(end,:)';
    dPhidx = (Phi_pert-Phi0)/norm(theta(j),2);
    Jf_FD(:,j) = dPhidx;
end
end

function [bodies] = nbody_init(labels)
%NBODY_INIT Initialize planetary data for n-body propagation
%   Given a set of labels of planets and/or barycentres, returns a
%   cell array populated with structures containing the body label and the
%   associated gravitational constant.
%
%
% Author
%   Name: ALESSANDRO 
%   Surname: MORSELLI
%   Research group: DART
%   Department: DAER
%   University: Politecnico di Milano 
%   Creation: 26/09/2021
%   Contact: alessandro.morselli@polimi.it
%   Copyright: (c) 2021 A. Morselli, Politecnico di Milano. 
%                  All rights reserved.
%
%
% Notes:
%   This material was prepared to support the course 'Satellite Guidance
%   and Navigation', AY 2021/2022.
%
%
% Inputs:
%   labels : [1,n] cell-array with object labels
%
% Outputs:
%   bodies : [1,n] cell-array with struct elements containing the following
%                  fields
%                  |
%                  |--bodies{i}.name -> body label
%                  |--bodies{i}.GM   -> gravitational constant [km**3/s**2]
%
%
% Prerequisites:
%   - MICE (Matlab SPICE)
%   - Populated kernel pool (PCK kernels)
%

% Initialize output
bodies = cell(size(labels));

% Loop over labels
for i = 1:length(labels)
    % Store body label
    bodies{i}.name = labels{i};
    % Store body gravitational constant
    bodies{i}.GM   = cspice_bodvrd(labels{i}, 'GM', 1);
end

end

function dy = ode_orbit(~, y, mu)
%
%     ode_orbit - 2-body problem system to be integrated
%
%     DESCRIPTION:
%       Function to be provided to the ODE integrator for the integration
%       of an orbit
%
%     PROTOTYPE:
%        dy = ode_orbit(~, y, mu)
%          
%     INPUT:
%       t              Integration time
%       y              Integration vector (the one with initial conditions)
%       mu             Gravitational constant
%     
%     OUTPUT:
%       dy             Vector to be integrated
%     
%     CALLED FUNCTIONS:
%       -
% 
%     LAST UPDATED:
%      6/11/2021
%
%     CREATED BY:
%      Bardazzi Niccolò

r = y(1:3); v = y(4:6);
dy(1:3,1) = v;           
dy(4:6,1) = -mu*r/norm(r,2).^3;  
end

function dy = ode_orbit_target(~, y, data)
%
%     ode_orbit_target - 2-body problem system to be integrated
%
%     DESCRIPTION:
%       Function to be provided to the ODE integrator for the integration
%       of an orbit
%
%     PROTOTYPE:
%        dy = ode_orbit_target(~, y, data)
%          
%     INPUT:
%       t              Integration time
%       y              Integration vector (the one with initial conditions)
%       data           Struct with fields:
%                      - data.TT: scaling in time
%                      - body{}: struct with the 2 planets name, example:
%                        {'Earth'}, {'Mars'}
%                      - data.LL: scaling in position
%                      - data.frame: 'J2000' 
%                      - data.mu: gravitational constant of the main
%                        attractor
%                      - data.turn: option in plotting trajectories
%                      - data.x_target: target position, useful as the
%                        event of the integration
%     
%     OUTPUT:
%       dy             Vector to be integrated
%     
%     CALLED FUNCTIONS:
%       -
%      useful with the event function: target_reached
% 
%     LAST UPDATED:
%      6/11/2021
%
%     CREATED BY:
%      Bardazzi Niccolò

mu = data.mu;
r = y(1:3); v = y(4:6);
dy(1:3,1) = v;           
dy(4:6,1) = -mu*r/norm(r,2).^3;  
end

function [value, isterminal, direction] = target_reached(~, y, data)
%
%     target_reached - Event of integration
%
%     DESCRIPTION:
%       Event of the integration when a certain point has been reached
%
%     PROTOTYPE:
%        [value, isterminal, direction] = target_reached(~, y, data)
%          
%     INPUT:
%       t              Initial state: [correct position; guess in velocity] 
%       y              Integration vector (the one with initial conditions)
%       data           Struct with fields:
%                      - data.TT: scaling in time
%                      - body{}: struct with the 2 planets name, example:
%                        {'Earth'}, {'Mars'}
%                      - data.LL: scaling in position
%                      - data.frame: 'J2000' 
%                      - data.mu: gravitational constant of the main
%                        attractor
%                      - data.turn: option in plotting trajectories
%                      - data.x_target: target position, useful as the
%                        event of the integration
%     
%     OUTPUT:
%       value          Function whose zero has been searched
%       isterminal     = 1, stop at the event
%       direction      = 0, generical
%     
%     CALLED FUNCTIONS:
%       -
% 
%     LAST UPDATED:
%      6/11/2021
%
%     CREATED BY:
%      Bardazzi Niccolò

    x_target = data.x_target;
    value = norm(x_target-y(1:3));
    isterminal = 1;
    direction = 0;    
end

%% Lambert's solver
function [A,P,E,ERROR,VI,VF,TPAR,THETA] = lambertMR(RI,RF,TOF,MU,orbitType,Nrev,Ncase,optionsLMR)

% lambertMR.m - Lambert's problem solver for all possible transfers
%   (multi-revolution transfer included).
%
% PROTOTYPE:
%   [A,P,E,ERROR,VI,VF,TPAR,THETA] = lambertMR(RI,RF,TOF,MU,orbitType,Nrev,Ncase,optionsLMR)
%
% DESCRIPTION:
%   Lambert's problem solver for all possible transfers:
%       1- zero-revolution (for all possible types of orbits: circles, ellipses,
%       	parabolas and hyperbolas)
%       2- multirevolution case
%       3- inversion of the motion
%
%   1- ZERO-REVOLUTION LAMBERT'S PROBLEM
%
%   For the solution of Lambert's problem with number of revolution = 0 the
%   subroutine by Chris D'Souza is included here.
%   This subroutine is a Lambert algorithm which given two radius vectors
%   and the time to get from one to the other, it finds the orbit
%   connecting the two. It solves the problem using a new algorithm
%   developed by R. Battin. It solves the Lambert problem for all possible
%   types of orbits (circles, ellipses, parabolas and hyperbolas).
%   The only singularity is for the case of a transfer angle of 360 degrees,
%   which is a rather obscure case.
%   It computes the velocity vectors corresponding to the given radius
%   vectors except for the case when the transfer angle is 180 degrees
%   in which case the orbit plane is ambiguous (an infinite number of
%   transfer orbits exist).
% 
%   2- MULTIREVOLUTION LAMBERT'S PROBLEM
%
%   For the solution of Lambert's problem with Nrev>0 number of revolution,
%   Battin's formulation has been extended to accomodate N-revolution
%   transfer orbits, by following the paper: "Using Battin Mathod to obtain 
%   Multiple-revolution Lambert's Solutions" by Shen and Tsiotras.
%
%   When Nrev>0 the possible orbits are just ellipses.
%   If 0<=Nrev<=Nmax, there are two Nrev-revolution transfer orbits.
%   These two transfer orbits have different semi-major axis and they may 
%   be all combinations of large-e and small-e transfer orbits.
%   The Original Successive Substitution Method by Battin converges to one
%   of the two possible solution with a viable initial guest, however it
%   diverges from the other one. Then a Reversed Successive Substitution is
%   used to converge to the second solution.
%   A procedure is implemented in order to guarantee to provide initial
%   guesses in the convergence region. If Nrew exceeds the maximum number
%   of revolution an ERROR is given:
%   warning('off','lambertMR:SuccessiveSubstitutionDiverged') to take out
%   the warnings or use optionsLMR(1) = 0.
% 
%   3- INVERSION OF THE MOTION
% 
%   Direct or retrograde option can be selected for the transfer.
%   
%   The algorithm computes the semi-major axis, the parameter (semi-latus 
%   rectum), the eccentricity and the velocity vectors.
% 
%   NOTE: If ERROR occurs or the 360 or 180 degree transfer case is 
%   encountered. 
%
% INPUT:
%	RI[3]           Vector containing the initial position in Cartesian
%                   coordinates [L].
%	RF[3]           Vector containing the final position vector in
%                   Cartesian coordinates [L].
%	TOF[1]          Transfer time, time of flight [T].
%  	MU[1]           Planetary constant of the planet (mu = mass * G) [L^3/T^2]
%	orbitType[1]    Logical variable defining whether transfer is
%                       0: direct transfer from R1 to R2 (counterclockwise)
%                       1: retrograde transfer from R1 to R2 (clockwise)
%	Nrev[1]         Number of revolutions.
%                   if Nrev = 0 ZERO-REVOLUTION transfer is calculated
%                   if Nrev > 0 two transfers are possible. Ncase should be
%                          defined to select one of the two.
%	Ncase[1]        Logical variable defining the small-a or large-a option
%                   in case of Nrev>0:
%                       0: small-a option
%                       1: large-a option
%	optionsLMR[1]	lambertMR options:
%                    optionsLMR(1) = display options:
%                                    0: no display
%                                    1: warnings are displayed only when
%                                       the algorithm does not converge
%                                    2: full warnings displayed
%
% OUTPUT:
%	A[1]        Semi-major axis of the transfer orbit [L].
% 	P[1]        Semi-latus rectum of the transfer orbit [L].
%  	E[1]        Eccentricity of the transfer orbit.
%	ERROR[1]	Error flag
%                   0:	No error
%                   1:	Error, routine failed to converge
%                   -1:	180 degrees transfer
%                   2:  360 degrees transfer
%                   3:  the algorithm doesn't converge because the number 
%                       of revolutions is bigger than Nrevmax for that TOF
%                   4:  Routine failed to converge, maximum number of
%                       iterations exceeded.
%	VI[3]       Vector containing the initial velocity vector in Cartesian
%               coordinates [L/T].
%	VT[1]		Vector containing the final velocity vector in Cartesian
%               coordinates [L/T].
%	TPAR[1] 	Parabolic flight time between RI and RF [T].
%	THETA[1]	Transfer angle [radians].
%
% NOTE: The semi-major axis, positions, times, and gravitational parameter
%       must be in compatible units.
%
% CALLED FUNCTIONS:
%   qck, h_E (added at the bottom of this file)
%
% REFERENCES:
%   - Shen and Tsiotras, "Using Battin method to obtain Multiple-Revolution
%       Lambert's solutions".
%   - Battin R., "An Introduction to the Mathematics and Methods of
%       Astrodynamics, Revised Edition", 1999.
%
% FUTURE DEVELOPMENT:
%   - 180 degrees transfer indetermination
%   - 360 degrees transfer singularity
%   - Nmax number of max revolution for a given TOF:
%     work in progress - Camilla Colombo
%
% ORIGINAL VERSION:
%   Chris D'Souza, 20/01/1989, MATLAB, lambert.m
%       verified by Darrel Monroe, 10/25/90
%       - Labert.m solved only direct transfer, without multi-revolution
%         option
%
% AUTHOR:
%   Camilla Colombo, 10/11/2006, MATLAB, lambertMR.m
%
% CHANGELOG:
%   13/11/2006, Camilla Colombo: added ERROR = 3 if Nrev > NrevMAX
%	21/11/2006, Camilla Colombo: added another case of ERROR = 3 (index
%   	N3) corresponding to the limit case when small-a solution = large-a
%       solution. No solution is given in this case.
%	06/08/2007, Camilla Colombo: optionsLMR added as an input
%	28/11/2007, Camilla Colombo: minor changes
%   29/01/2009, Matteo Ceriotti:
%       - Introduced variable for maximum number of iterations nitermax.
%       - Corrected final check on maximum number of iterations exceeded, from
%           "==" to ">=" (if N1 >= nitermax || N >= nitermax).
%       - Increased maxumum number of iterations to 2000, not to lose some
%           solutions.
%       - In OSS loop, added check for maximum number of iterations exceeded,
%           which then sets checkNconvOSS = 0.
%       - Changed the way of coumputing X given Y1 in RSS. Now the
%           Newton-Raphson method with initial guess suggested by Shen,
%           Tsiotras is used. This should guarantee convergence without the
%           need of an external zero finder (fsolve).
%       - Changed absolute tolerance into relative tolerance in all loops X0-X.
%           Now the condition is: while "abs(X0-X) >= abs(X)*TOL+TOL".
%       - Added return immediately when any error is detected.
%       - Moved check on 4*TOF*LAMBDA==0 after computing LAMBDA.
%       - Moved check on THETA==0 || THETA==2*PI after computing THETA.
%       - Added error code 4 (number of iterations exceeded).
%       - Removed variable Nwhile, as not strictly needed.
%       - Removed variable PIE=pi.
%   29/01/2009, REVISION: Matteo Ceriotti
%   21/07/2009, Matteo Ceriotti, Camilla Colombo:
%       added condition to detect case 180 degrees transfer indetermination
%   30/01/2010, Camilla Colombo: Header and function name in accordance
%       with guidlines.
%
% Note: Please if you have got any changes that you would like to be done,
%   do not change the function, please contact the author.
%
% -------------------------------------------------------------------------

% Check inputs
if nargin < 8
    optionsLMR = 0;
    if nargin < 6
        Nrev = 0;
        if nargin < 5
            orbitType = 0;
            if nargin < 4
                error('Not enough input arguments. See lambertMR.');
            end
        end
    end
end

nitermax = 2000; % Maximum number of iterations for loops
TOL = 1e-14;

TWOPI=2*pi;

% Reset
A=0;P=0;E=0;VI=[0,0,0];VF=[0,0,0];

% ----------------------------------
% Compute the vector magnitudes and various cross and dot products

RIM2   = dot(RI,RI);
RIM    = sqrt(RIM2);
RFM2   = dot(RF,RF);
RFM    = sqrt(RFM2);
CTH    = dot(RI,RF)/(RIM*RFM);
CR     = cross(RI,RF);
STH    = norm(CR)/(RIM*RFM);

% Choose angle for up angular momentum
switch orbitType
    case 0 % direct transfer
        if CR(3) < 0 
            STH = -STH;
        end
    case 1 % retrograde transfer
        if CR(3) > 0 
            STH = -STH;
        end
    otherwise
		error('%d is not an allowed orbitType',orbitType);
end
        
THETA  = qck(atan2(STH,CTH));
% if abs(THETA - pi) >= 0.01
if THETA == TWOPI || THETA==0
    ERROR = 2;
    A=0; P=0; E=0; VI=[0,0,0]; VF=[0,0,0]; TPAR=0; THETA=0;
    return
end

B1     = sign(STH); if STH == 0; B1 = 1; end;

% ----------------------------------
% Compute the chord and the semi-perimeter

C= sqrt(RIM2 + RFM2 - 2*RIM*RFM*CTH);
S= (RIM + RFM + C)/2;
BETA   = 2*asin(sqrt((S-C)/S));
PMIN   = TWOPI*sqrt(S^3/(8*MU));
TMIN   = PMIN*(pi-BETA+sin(BETA))/(TWOPI);
LAMBDA = B1*sqrt((S-C)/S);

if 4*TOF*LAMBDA == 0 || abs((S-C)/S) < TOL
    ERROR = -1;
    A=0; P=0; E=0; VI=[0,0,0]; VF=[0,0,0]; TPAR=0; THETA=0;
    return
end

% ----------------------------------
% Compute L carefully for transfer angles less than 5 degrees

if THETA*180/pi <= 5
   W   = atan((RFM/RIM)^.25) - pi/4;
   R1  = (sin(THETA/4))^2;
   S1  = (tan(2*W))^2;
   L   = (R1+S1)/(R1+S1+cos(THETA/2));
else
   L   = ((1-LAMBDA)/(1+LAMBDA))^2;
end

M= 8*MU*TOF^2/(S^3*(1+LAMBDA)^6);
TPAR   = (sqrt(2/MU)/3)*(S^1.5-B1*(S-C)^1.5);
L1     = (1 - L)/2;

CHECKFEAS = 0;
N1 = 0;
N = 0;

if Nrev == 0
    % ----------------------------------
    % Initialize values of y, n, and x

    Y= 1;
    N= 0;
    N1=0;
    ERROR  = 0;
    % CHECKFEAS=0;

    if (TOF-TPAR) <= 1e-3
        X0  = 0;
    else
        X0  = L;
    end

    X= -1.e8;

    % ----------------------------------
    % Begin iteration
    
    % ---> CL: 26/01/2009, Matteo Ceriotti: 
    %       Changed absolute tolerance into relative tolerance here below.
    while (abs(X0-X) >= abs(X)*TOL+TOL) && (N <= nitermax)
        N   = N+1;
        X   = X0;
        ETA = X/(sqrt(1+X) + 1)^2;
        CHECKFEAS=1;

        % ----------------------------------
        % Compute x by means of an algorithm devised by
        % Gauticci for evaluating continued fractions by the
        % 'Top Down' method
        
        DELTA = 1;
        U     = 1;
        SIGMA = 1;
        M1    = 0;

        while abs(U) > TOL && M1 <= nitermax
            M1    = M1+1;
            GAMMA = (M1 + 3)^2/(4*(M1+3)^2 - 1);
            DELTA = 1/(1 + GAMMA*ETA*DELTA);
            U     = U*(DELTA - 1);
            SIGMA = SIGMA + U;
        end

        C1 = 8*(sqrt(1+X)+1)/(3+1/(5 + ETA + (9*ETA/7)*SIGMA));

        % ----------------------------------
        % Compute H1 and H2
        
        if N == 1
            DENOM = (1 + 2*X + L)*(3*C1 + X*C1 +4*X);
            H1 = (L+X)^2*(C1 + 1 + 3*X)/DENOM;
            H2 = M*(C1+X-L)/DENOM;
        else
            QR = sqrt(L1^2 + M/Y^2);
            XPLL = QR - L1;
            LP2XP1 = 2*QR;
            DENOM = LP2XP1*(3*C1 + X*C1+4*X);
            H1 = ((XPLL^2)*(C1 + 1 + 3*X))/DENOM;
            H2 = M*(C1+X-L)/DENOM;
        end
        
        B = 27*H2/(4*(1+H1)^3);
        U = -B/(2*(sqrt(B+1)+1));

        % ----------------------------------
        % Compute the continued fraction expansion K(u)
        % by means of the 'Top Down' method
        
        % Y can be computed finding the roots of the formula and selecting
        % the real one:
        % y^3 - (1+h1)*y^2 - h2 = 0     (7.113) Battin
        %
        % Ycami_ = roots([1 -1-H1 0 -H2])
        % kcami = find( abs(imag(Ycami_)) < eps );
        % Ycami = Ycami_(kcami)

        DELTA = 1;
        U0 = 1;
        SIGMA = 1;
        N1 = 0;

        while N1 < nitermax && abs(U0) >= TOL
            if N1 == 0
                GAMMA = 4/27;
                DELTA = 1/(1-GAMMA*U*DELTA);
                U0 = U0*(DELTA - 1);
                SIGMA = SIGMA + U0;
            else
                for I8 = 1:2
                    if I8 == 1
                        GAMMA = 2*(3*N1+1)*(6*N1-1)/(9*(4*N1 - 1)*(4*N1+1));
                    else
                        GAMMA = 2*(3*N1+2)*(6*N1+1)/(9*(4*N1 + 1)*(4*N1+3));
                    end
                    DELTA = 1/(1-GAMMA*U*DELTA);
                    U0 = U0*(DELTA-1);
                    SIGMA = SIGMA + U0;
                end
            end

            N1 = N1 + 1;
        end

        KU = (SIGMA/3)^2;
        Y = ((1+H1)/3)*(2+sqrt(B+1)/(1-2*U*KU));    % Y = Ycami
        
        X0 = sqrt(((1-L)/2)^2+M/Y^2)-(1+L)/2;
        % fprintf('n= %d, x0=%.14f\n',N,X0);
    end
    
% MULTIREVOLUTION
elseif (Nrev > 0) && (4*TOF*LAMBDA~=0) %(abs(THETA)-pi > 0.5*pi/180)

    checkNconvRSS = 1;
    checkNconvOSS = 1;
    N3 = 1;
    
    while N3 < 3
        
        if Ncase == 0 || checkNconvRSS == 0

            % - Original Successive Substitution -
            % always converges to xL - small a

            % ----------------------------------
            % Initialize values of y, n, and x
            
            Y= 1;
            N= 0;
            N1=0;
            ERROR = 0;
            % CHECKFEAS = 0;
%             if (TOF-TPAR) <= 1e-3
%                 X0 = 0;
%             else
            if checkNconvOSS == 0
                X0 = 2*X0;
                checkNconvOSS = 1;
                % see p. 11 USING BATTIN METHOD TO OBTAIN 
                % MULTIPLE-REVOLUTION LAMBERT'S SOLUTIONS - Shen, Tsiotras
            elseif checkNconvRSS == 0;
                % X0 is taken from the RSS
            else
                X0 = L;
            end

            X = -1.e8;

            % ----------------------------------
            % Begin iteration
            			
            % ---> CL: 26/01/2009,Matteo Ceriotti 
            %   Changed absolute tolerance into relative tolerance here
            %   below.
            while (abs(X0-X) >= abs(X)*TOL+TOL) && (N <= nitermax)
                N   = N+1;
                X   = X0;
                ETA = X/(sqrt(1+X) + 1)^2;
                CHECKFEAS = 1;

                % ----------------------------------
                % Compute x by means of an algorithm devised by
                % Gauticci for evaluating continued fractions by the
                % 'Top Down' method
                

                DELTA = 1;
                U     = 1;
                SIGMA = 1;
                M1    = 0;

                while abs(U) > TOL && M1 <= nitermax
                    M1    = M1+1;
                    GAMMA = (M1 + 3)^2/(4*(M1+3)^2 - 1);
                    DELTA = 1/(1 + GAMMA*ETA*DELTA);
                    U     = U*(DELTA - 1);
                    SIGMA = SIGMA + U;
                end

                C1 = 8*(sqrt(1+X)+1)/(3+1/(5 + ETA + (9*ETA/7)*SIGMA));

                % ----------------------------------
                % Compute H1 and H2
                
                if N == 1
                    DENOM = (1 + 2*X + L)*(3*C1 + X*C1 +4*X);
                    H1 = (L+X)^2*(C1 + 1 + 3*X)/DENOM;
                    H2 = M*(C1+X-L)/DENOM;
                else
                    QR = sqrt(L1^2 + M/Y^2);
                    XPLL = QR - L1;
                    LP2XP1 = 2*QR;
                    DENOM = LP2XP1*(3*C1 + X*C1+4*X);
                    H1 = ((XPLL^2)*(C1 + 1 + 3*X))/DENOM;
                    H2 = M*(C1+X-L)/DENOM;
                end

                H3 = M*Nrev*pi/(4*X*sqrt(X));
                H2 = H3+H2;

                B = 27*H2/(4*(1+H1)^3);
                U = -B/(2*(sqrt(B+1)+1));

                % ----------------------------------
                % Compute the continued fraction expansion K(u)
                % by means of the 'Top Down' method
                
                % Y can be computed finding the roots of the formula and selecting
                % the real one:
                % y^3 - (1+h1)*y^2 - h2 = 0     (7.113) Battin
                %
                % Ycami_ = roots([1 -1-H1 0 -H2])
                % kcami = find( abs(imag(Ycami_)) < eps );
                % Ycami = Ycami_(kcami)

                DELTA = 1;
                U0 = 1;
                SIGMA = 1;
                N1 = 0;

                while N1 < nitermax && abs(U0) >= TOL
                    if N1 == 0
                        GAMMA = 4/27;
                        DELTA = 1/(1-GAMMA*U*DELTA);
                        U0 = U0*(DELTA - 1);
                        SIGMA = SIGMA + U0;
                    else
                        for I8 = 1:2
                            if I8 == 1
                                GAMMA = 2*(3*N1+1)*(6*N1-1)/(9*(4*N1 - 1)*(4*N1+1));
                            else
                                GAMMA = 2*(3*N1+2)*(6*N1+1)/(9*(4*N1 + 1)*(4*N1+3));
                            end
                            DELTA = 1/(1-GAMMA*U*DELTA);
                            U0 = U0*(DELTA-1);
                            SIGMA = SIGMA + U0;
                        end
                    end

                    N1 = N1 + 1;
                end

                KU = (SIGMA/3)^2;
                Y = ((1+H1)/3)*(2+sqrt(B+1)/(1-2*U*KU));	% Y = Ycami
                if Y > sqrt(M/L)
                    if optionsLMR(1) == 2
                        warning('lambertMR:SuccessiveSubstitutionDiverged',...
                                ['Original Successive Substitution is diverging\n'...
                                '-> Reverse Successive Substitution used to find the proper XO.\n']);
                    end
                    checkNconvOSS = 0;
                    break
                end
                
                X0 = sqrt(((1-L)/2)^2+M/Y^2)-(1+L)/2;
                % fprintf('N: %d X0: %.14f\n',N,X0);
            end
            
            % When 2 solutions exist (small and big a), the previous loop
            % must either converge or diverge because Y > sqrt(M/L) at some
            % point. Thus, the upper bound on the number of iterations
            % should not be necessary. Though, nothing can be said in the
            % case tof<tofmin and so no solution exist. In this case, an
            % upper bound on number of iterations could be needed.
            
            if N >= nitermax % Checks if previous loop ended due to maximum number of iterations
                if optionsLMR(1) == 2
                    warning('lambertMR:SuccessiveSubstitutionExceedMaxIter',...
                            ['Original Successive Substitution exceeded max number of iteration\n'...
                            '-> Reverse Successive Substitution used to find the proper XO.\n']);
                end
                checkNconvOSS = 0;
            end
        end
        if (Ncase == 1 || checkNconvOSS == 0) && ~(checkNconvRSS == 0 && checkNconvOSS == 0)

            % - Reverse Successive Substitution -
            % always converges to xR - large a

            % ----------------------------------
            % Initialize values of y, n, and x
            
            N = 0;
            N1 = 0;
            ERROR  = 0;
            % CHECKFEAS=0;
            if checkNconvRSS == 0;
                X0 = X0/2; % XL/2
                checkNconvRSS = 1;
                % see p. 11 USING BATTIN METHOD TO OBTAIN 
                % MULTIPLE-REVOLUTION LAMBERT'S SOLUTIONS - Shen, Tsiotras
            elseif checkNconvOSS == 0
                % X0 is taken from the OSS
            else
                X0 = L;
            end

            X = -1.e8;

            % ----------------------------------
            % Begin iteration
            
            % ---> CL: 26/01/2009, Matteo Ceriotti
            %   Changed absolute tolerance into relative tolerance here
            %   below.
            while (abs(X0-X) >= abs(X)*TOL+TOL) && (N <= nitermax)
                N = N+1;
                X = X0;
                CHECKFEAS=1;

                Y = sqrt(M/((L+X)*(1+X))); % y1 in eq. (8a) in Shen, Tsiotras

                if Y < 1
                    if optionsLMR(1) == 2
                        warning('lambertMR:SuccessiveSubstitutionDiverged',...
                                ['Reverse Successive Substitution is diverging\n' ...
                                '-> Original Successive Substitution used to find the proper XO.\n']);
                    end
                    checkNconvRSS = 0;
                    break
                end
                
                % ---> CL: 27/01/2009, Matteo Ceriotti
                %   This is the Newton-Raphson method suggested by USING
                %   BATTIN METHOD TO OBTAIN MULTIPLE-REVOLUTION LAMBERT'S
                %   SOLUTIONS - Shen, Tsiotras
                
                % To assure the Newton-Raphson method to be convergent
                Erss = 2*atan(sqrt(X));
                while h_E(Erss,Y,M,Nrev) < 0
                    Erss = Erss/2;
                end
                
                Nnew = 1;
                Erss_old = -1.e8;
                
                % The following Newton-Raphson method should always
                % converge, given the previous first guess choice,
                % according to the paper. Therefore, the condition on
                % number of iterations should not be neccesary. It could be
                % necessary for the case tof < tofmin.
                while (abs(Erss-Erss_old) >= abs(Erss)*TOL+TOL) && Nnew < nitermax
                    Nnew = Nnew+1;
                    [h, dh] = h_E(Erss,Y,M,Nrev);
                    Erss_old = Erss;
                    Erss = Erss - h/dh;
                    % fprintf('Nnew: %d Erss: %.16f h_E: %.16f\n',Nnew,Erss,h);
                end
                if Nnew >= nitermax
                    if optionsLMR(1) ~= 0
                        warning('lambertMR:NewtonRaphsonIterExceeded', 'Newton-Raphson exceeded max iterations.\n');
                    end
                end
                X0 = tan(Erss/2)^2;
            end
        end
        if checkNconvOSS == 1 && checkNconvRSS == 1
            break
        end
        
        if checkNconvRSS == 0 && checkNconvOSS == 0
            if optionsLMR ~=0
                warning('lambertMR:SuccessiveSubstitutionDiverged',...
                        ['Both Original Successive Substitution and Reverse ' ...
                        'Successive Substitution diverge because Nrev > NrevMAX.\n' ...
                        'Work in progress to calculate NrevMAX.\n']);
            end
            ERROR = 3;
            A=0; P=0; E=0; VI=[0,0,0]; VF=[0,0,0]; TPAR=0; THETA=0;
            return
        end
        
        N3 = N3+1;
    end
    
    if N3 == 3
        if optionsLMR ~=0
            warning('lambertMR:SuccessiveSubstitutionDiverged',...
                    ['Either Original Successive Substitution or Reverse ' ...
                    'Successive Substitution is always diverging\n' ...
                    'because Nrev > NrevMAX or because large-a solution = small-a solution (limit case).\n' ...
                    'Work in progress to calculate NrevMAX.\n']);
        end
        ERROR = 3;
        A=0; P=0; E=0; VI=[0,0,0]; VF=[0,0,0]; TPAR=0; THETA=0;
        return
    end
end

% ----------------------------------
% Compute the velocity vectors

if CHECKFEAS == 0
    ERROR = 1;
    A=0; P=0; E=0; VI=[0,0,0]; VF=[0,0,0]; TPAR=0; THETA=0;
    return
end

if N1 >= nitermax || N >= nitermax
    ERROR = 4;
    if optionsLMR ~=0
        disp('Lambert algorithm has not converged, maximum number of iterations exceeded.');
    end
    A=0; P=0; E=0; VI=[0,0,0]; VF=[0,0,0]; TPAR=0; THETA=0;
    return
end

CONST = M*S*(1+LAMBDA)^2;
A = CONST/(8*X0*Y^2);

R11 = (1 + LAMBDA)^2/(4*TOF*LAMBDA);
S11 = Y*(1 + X0);
T11 = (M*S*(1+LAMBDA)^2)/S11;

VI(1:3) = -R11*(S11*(RI(1:3)-RF(1:3))-T11*RI(1:3)/RIM);
VF(1:3) = -R11*(S11*(RI(1:3)-RF(1:3))+T11*RF(1:3)/RFM);

P = (2*RIM*RFM*Y^2*(1+X0)^2*sin(THETA/2)^2)/CONST;
E = sqrt(1 - P/A);

return

% -------------------------------------------------------------------------



function [angle] = qck(angle)

% qck.m - Reduce an angle between 0 and 2*pi
%
% PROTOTYPE:
%   [angle]=qck(angle)
%
% DESCRIPTION:
%   This function takes any angle and reduces it, if necessary,
% 	so that it lies in the range from 0 to 2 PI radians.
% 
% INPUTS:
%   ANGLE[1]    Angle to be reduced (in radians)
% 
% OUTPUTS:
%   QCK[1]      The angle reduced, if necessary, to the range
%               from 0 to 2 PI radians (in radians)
% 
% CALLED FUNCTIONS:
%   pi (from MATLAB)
%
% AUTHOR:
%   W.T. Fowler, July, 1978
%
% CHANGELOG:
%   8/20/90, REVISION: Darrel Monroe
%
% -------------------------------------------------------------------------

twopi = 2*pi;
 
diff = twopi * (fix(angle/twopi) + min([0,sign(angle)]));

angle = angle -diff;

return

end

% -------------------------------------------------------------------------

function [h, dh] = h_E(E, y, m, Nrev)

% h_E.m - Equation of multirevolution Lambert's problem h = h(E).
%
% PROTOTYPE:
%   [h, dh] = h_E(E, y, m, Nrev)
%
% DESCRIPTION:
%   Equation of multirevolution Lambert's problem:
%   h(E) = (Nrev*pi + E - sin(E)) / tan(E/2)^3 - 4/m * (y^3 - y^2)
%   See: "USING BATTIN METHOD TO OBTAIN MULTIPLE-REVOLUTION LAMBERT'S 
%      SOLUTIONS", Shen, Tsiotras, pag. 12
%
% INPUT
%   E, y, m, Nrev   See paper for detailed description.
%
% OUTPUT
%   h               Value of h(E).
%   dh              Value of dh(E)/dE.
%
% ORIGINAL VERSION:
%   Camilla Colombo, 20/02/2006, MATLAB, cubicN.m
%
% AUTHOR:
%   Matteo Ceriotti, 27/01/2009
%   - changed name of cubicN.m and added at the bottom of lambertMR.m file
%
% -------------------------------------------------------------------------

tanE2 = tan(E/2);
h = (Nrev*pi + E - sin(E)) / tanE2^3 - 4/m * (y^3 - y^2);

if nargout > 1  % two output arguments
    % h'(E)
    dh = (1-cos(E))/tanE2^3 - 3/2*(Nrev*pi+E-sin(E))*sec(E/2)^2 / tanE2^4;
end

return
end


end
% Spacecraft Guidance & Control (2021/2022)
% Assignment # 2
% Author: Niccolò Bardazzi 10800456

% Loading kernels
cspice_furnsh('assignment02.tm');
addpath(genpath(pwd))

%% Ex 1
clearvars; close all; clc; 

mu = cspice_bodvrd('EARTH','GM', 3);

spacecraft.ID = 28922;
spacecraft.Name = 'GIOVE-A';

x_mean_0 = [ 23721.610
             17903.673
              -49.918
             -1.150987
              1.529718
              3.122389 ];

P_0 = [  2.6e-2 1.4e-2 -1.8e-3 0 0 0
         1.4e-2 1.8e-2  2.3e-3 0 0 0
        -1.8e-3 2.3e-3  1.0e-2 0 0 0
        zeros(3,3) eye(3)*1.6e-7];

et0 = cspice_str2et('2021-11-19T14:25:39.652');
etf = cspice_str2et('2021-11-23T14:00:00');

options_ode = odeset('reltol', 1e-12, 'abstol', 1e-12);
N = floor((etf-et0)/60);
t = zeros(1,N+1);
t(1) = et0;
et = cspice_str2et('2021-11-19T14:26:00');
for i = 2:N+1
    t(i) = et;
    et = et+60;
end

% t = linspace(et0, etf, N);
[~, x_ref] = ode113(@(t,x) ode_orbit(t,x,mu), t, x_mean_0, options_ode); 

% Stations places
station(1).Name = 'Milano';
station(2).Name = 'Wellington';
station(3).Name = 'La-Silla';

% Latitude
station(1).LAT = 45.50122;
station(2).LAT = -41.28437;
station(3).LAT = -29.26117;

% Longitude
station(1).LON = 9.15461;
station(2).LON = 174.76697;
station(3).LON = -70.73133;

% Altitude
station(1).ALT = 20; % m
station(2).ALT = 117; % m
station(3).ALT = 2400; % m

% Type of sensor
station(1).type = 'Radar'; % monostatic
station(2).type = 'Optical';
station(3).type = 'Optical';

% Noise
station(1).noise = 0.1;          station(1).noise_range = 0.01; % km
station(2).noise = 0.0005;       station(2).noise_range = 0.01; % km
station(3).noise = 0.001;        station(3).noise_range = 0.01; % km

% FOV
station(1).FOV = 6; 
station(2).FOV = 2; 
station(3).FOV = 1; 

% Minimum Elevation
station(1).min_El = 30; 
station(2).min_El = 20; 
station(3).min_El = 10; 


for i = 1:length(station)
    [s3a_azimuth, s3a_elevation, s3a_range, s3a_range_rate] = ...
        antenna_pointing(station(i).Name, t, x_ref');
    station(i).visibility = s3a_elevation >= deg2rad(station(i).min_El);

    num = 1;
    long = 1;
    for a = 1:length(station(i).visibility)
        if station(i).visibility(a) ~= 0
            station(i).windows{num,long} = t(a);
            long = long+1;
        else
            if long ~= 1
                num = num+1;
                long = 1;
            end
        end
    end

end

tformat = 'YYYY-MON-DD-HR:MN:SC.####::UTC'; 
last_epoch = zeros(3,1);
for i = 1:length(station)
    last_window = zeros(1,2);
    b = 1;
    while isempty(station(i).windows{end,b}) == 0
        last_window(b) = station(i).windows{end,b};
        b = b+1;
    end
    last_window_start = cspice_timout(last_window(1),tformat); 
    last_window_end = cspice_timout(last_window(end),tformat); 
    last_epoch(i) = last_window(end);
    fprintf(['\n Last window of ', station(i).Name, ' is from ',last_window_start,' to ', last_window_end,' \n \n']);
end


%% Linearized and unscented an Montecarlo
% Linearized approach
for i = 1:length(station)
    N = floor((last_epoch(i)-et0)/60+1);
    t_LinCov = linspace(et0, last_epoch(i), N);
    I = eye(6);
    yPhi0 = [x_mean_0; I(:)];
    [~, LinCov_prop] = ode113(@(t,x) ode_STM(t,x, mu), t_LinCov, yPhi0, options_ode);
    station(i).LinCov.mean_end = LinCov_prop(end,1:6)';
    STM = reshape(LinCov_prop(end, 7:42), [6, 6]);
    station(i).LinCov.cov = STM*P_0*STM';
end

% Unscented approach
alpha = 1.e-2;
beta = 2;
lambda = alpha^2*(length(P_0))-length(P_0);
W0_m = lambda/(length(P_0)+lambda);
Wi_m = 1/(2*(length(P_0)+lambda));
W0_c = lambda/(length(P_0)+lambda)+(1-alpha^2+beta);
Wi_c = Wi_m;
for i = 1:length(station)
    sigmas = sigma_points_weighted(P_0, x_mean_0, lambda); % Sigma points generation
    sz = size(sigmas);
    y = zeros(6,sz(2));
    N = floor((last_epoch(i)-et0)/60+1);
    t_UT = linspace(et0, last_epoch(i), N);
    sum1 = zeros(6,1);
    sum2 = zeros(6,6);
    for j = 1:sz(2)
        [~, UT_prop1] = ode113(@(t,x) ode_orbit(t,x,mu), t_UT, sigmas(:,j), options_ode); 
        y(:,j) = UT_prop1(end,:);
        % Mean  
        if j == 1
            sum1 = UT_prop1(end,:)'*W0_m+sum1;  
        else
            sum1 = UT_prop1(end,:)'*Wi_m+sum1;  
        end
    end
    station(i).UT.mean_end = sum1;
    for j = 1:sz(2)
        % Covariance
        if j == 1
            sum2 = (y(:,j)-station(i).UT.mean_end)*(y(:,j)-station(i).UT.mean_end)'*W0_c+sum2;  
        else
            sum2 = (y(:,j)-station(i).UT.mean_end)*(y(:,j)-station(i).UT.mean_end)'*Wi_c+sum2;  
        end
    end
    station(i).UT.cov = sum2;
end

% Montecarlo
n = 100; % # of points
R = mvnrnd(x_mean_0,P_0,n);
for i = 1:length(station)
    N = floor((last_epoch(i)-et0)/60+1);
    t_MC = linspace(et0, last_epoch(i), N);
    for j = 1:n
        [~, MC_prop] = ode113(@(t,x) ode_orbit(t,x,mu), t_MC, R(j,:), options_ode); 
        y(:,j) = MC_prop(end,:);
    end
    station(i).MC.samples = y;
    % Mean
    station(i).MC.mean_end = mean(y,2);
    % Covariance
    station(i).MC.cov = cov(y');
end

%% FOV percentage of samples with FOV pointing at reference trajectory
reci_station = zeros(3,3);
for i = 1:length(station)
    counter = 0;
    end_window = find(last_epoch(i) == t);  
    reci_station(:,i) = cspice_spkpos(station(i).Name, last_epoch(i), 'J2000', 'NONE', 'EARTH');
    for j = 1:n
        [~, ~, ~, ~, ECI2TOPO] = antenna_pointing(station(i).Name, last_epoch(i), station(i).MC.samples(:,j)); 
        rho = ECI2TOPO(1:3,1:3)*(station(i).MC.samples(1:3,j)- reci_station(:,i)); % in topocentric ref. frame
        rho_dir = rho/norm(rho);
        r_ref = x_ref(end_window,1:3);
        r_refTOPO = ECI2TOPO(1:3,1:3)*(r_ref'-reci_station(:,i));
        r_ref_dir = r_refTOPO/norm(r_refTOPO);
        alpha = acosd(dot(rho_dir, r_ref_dir));

        if alpha <= station(i).FOV/2
            counter = counter+1;
        end
    end
    station(i).percentage = counter/n*100;
end

%%
% pick the Station To Be Seen
stbs = 2; % 1: Milan
          % 2: Wellington
          % 3: La-Silla

fh = figure(4);
fh.WindowState = 'maximized';
% Mean values
scatter3(station(stbs).LinCov.mean_end(1), station(stbs).LinCov.mean_end(2), ...
    station(stbs).LinCov.mean_end(3),'LineWidth',1.3, ...
    'Marker','o' ,'Displayname', 'LinCov'), hold on
scatter3(station(stbs).UT.mean_end(1), station(stbs).UT.mean_end(2), ...
    station(stbs).UT.mean_end(3),'LineWidth',1.3, ...
    'Marker','o' ,'Displayname', 'Unscented Transform'), hold on
scatter3(station(stbs).MC.mean_end(1), station(stbs).MC.mean_end(2), ...
    station(stbs).MC.mean_end(3),'LineWidth',1.3, ...
    'Marker','o' ,'Displayname', 'Montecarlo'), hold on
% FOV visualization (cone)
R = linspace(0, norm(station(stbs).LinCov.mean_end(1:3))*sin(deg2rad(station(stbs).FOV)), 100);
[X, Y, Z] = cylinder2P(R, 100,reci_station(:,stbs)',station(stbs).LinCov.mean_end(1:3)');
surf(X,Y,Z,'LineStyle','none','FaceAlpha',0.3,'DisplayName',['FOV = ', num2str(station(stbs).FOV),'°']);
legend show, legend('AutoUpdate','off')
% Ellipsoids plot
h3 = plot_gaussian_ellipsoid(station(stbs).MC.mean_end(1:3), station(stbs).MC.cov(1:3,1:3), 3, 20);
set(h3,'facealpha',0.45);
set(h3,'LineStyle','none','FaceColor',[0.9290 0.6940 0.1250])
% Orbit
plot3(x_ref(:,1), x_ref(:,2), x_ref(:,3),'LineWidth',0.01,'Color','blue'), hold on
% Earth
RE = cspice_bodvrd('Earth','RADII',3);
[X, Y, Z] = sphere;
hSurface = surf(X*RE(1), Y*RE(2), Z*RE(3)); hold on
set(hSurface,'FaceColor',[0.5 0.5 0.5])
axis equal
grid on
title([spacecraft.Name, ' visibility from ', station(stbs).Name])
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
hold off
view(20,20)


fh = figure(5);
fh.WindowState = 'maximized';
% Mean values
scatter3(station(stbs).LinCov.mean_end(1), station(stbs).LinCov.mean_end(2), ...
    station(stbs).LinCov.mean_end(3),'LineWidth',1.3, ...
    'Marker','o' ,'Displayname', 'LinCov'), hold on
scatter3(station(stbs).UT.mean_end(1), station(stbs).UT.mean_end(2), ...
    station(stbs).UT.mean_end(3),'LineWidth',1.3, ...
    'Marker','o' ,'Displayname', 'Unscented Transform'), hold on
scatter3(station(stbs).MC.mean_end(1), station(stbs).MC.mean_end(2), ...
    station(stbs).MC.mean_end(3),'LineWidth',1.3, ...
    'Marker','o' ,'Displayname', 'Montecarlo'), hold on
% FOV visualization (cone)
R = linspace(0, norm(station(stbs).LinCov.mean_end(1:3)-reci_station(:,stbs))*sin(deg2rad(station(stbs).FOV)), 100);
[X, Y, Z] = cylinder2P(R, 100,reci_station(:,stbs)',station(stbs).LinCov.mean_end(1:3)');
surf(X,Y,Z,'LineStyle','none','FaceAlpha',0.3,'DisplayName',['FOV = ', num2str(station(stbs).FOV),'°']);
legend show, legend('AutoUpdate','off')
% Montecarlo points
max_dist = norm([X(end,end);Y(end,end);Z(end,end)]-station(stbs).LinCov.mean_end(1:3));
count1 = 0;
count2 = 0;
for j = 1:n
    if max_dist < norm(station(stbs).MC.samples(1:3,j)-station(stbs).LinCov.mean_end(1:3))
        color = 'red';
        if count1 == 0
            count1 = 1;
            legend('AutoUpdate','on')
            scatter3(station(stbs).MC.samples(1,j), station(stbs).MC.samples(2,j), station(stbs).MC.samples(3,j), 10, color, 'Marker','*','DisplayName','Outside FOV');
            legend('AutoUpdate','off')
        else
            scatter3(station(stbs).MC.samples(1,j), station(stbs).MC.samples(2,j), station(stbs).MC.samples(3,j), 10, color, 'Marker','*');
        end
    else
        color = [0.4940 0.1840 0.5560];
        if count2 == 0
            count2 = 1;
            legend('AutoUpdate','on')
            scatter3(station(stbs).MC.samples(1,j), station(stbs).MC.samples(2,j), station(stbs).MC.samples(3,j), 10, color, 'Marker','*','DisplayName','Inside FOV');
            legend('AutoUpdate','off')
        else
            scatter3(station(stbs).MC.samples(1,j), station(stbs).MC.samples(2,j), station(stbs).MC.samples(3,j), 10, color, 'Marker','*');
        end
    end
end
% Ellipsoids plot
h3 = plot_gaussian_ellipsoid(station(stbs).MC.mean_end(1:3), station(stbs).MC.cov(1:3,1:3), 3, 20);
set(h3,'facealpha',0.45);
set(h3,'LineStyle','none','FaceColor',[0.9290 0.6940 0.1250])
% Orbit
plot3(x_ref(:,1), x_ref(:,2), x_ref(:,3),'LineWidth',0.01,'Color','blue'), hold on
% Earth
RE = cspice_bodvrd('Earth','RADII',3);
[X, Y, Z] = sphere;
hSurface = surf(X*RE(1), Y*RE(2), Z*RE(3)); hold on
set(hSurface,'FaceColor',[0.5 0.5 0.5])
axis equal
grid on
title([spacecraft.Name, ' visibility from ', station(stbs).Name])
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
% Zooming from the previous figure
xlim([station(stbs).LinCov.mean_end(1)-sind(station(stbs).FOV)*norm(station(stbs).LinCov.mean_end), station(stbs).LinCov.mean_end(1)+sind(station(stbs).FOV)*norm(station(stbs).LinCov.mean_end)])
ylim([station(stbs).LinCov.mean_end(2)-sind(station(stbs).FOV)*norm(station(stbs).LinCov.mean_end), station(stbs).LinCov.mean_end(2)+sind(station(stbs).FOV)*norm(station(stbs).LinCov.mean_end)])
zlim([station(stbs).LinCov.mean_end(3)-.8*sind(station(stbs).FOV)*norm(station(stbs).LinCov.mean_end), station(stbs).LinCov.mean_end(3)+.8*sind(station(stbs).FOV)*norm(station(stbs).LinCov.mean_end)])
hold off

%% Ex 2
clearvars -except station last_epoch reci_station; close all; clc; 

mu = cspice_bodvrd('EARTH','GM', 3);

initial.mean = [ 23721.610
                 17903.673
                  -49.918
                 -1.150987
                  1.529718
                  3.122389 ];

initial.cov = [  2.6e-2 1.4e-2 -1.8e-3 0 0 0
                 1.4e-2 1.8e-2  2.3e-3 0 0 0
                -1.8e-3 2.3e-3  1.0e-2 0 0 0
                zeros(3,3) eye(3)*1.6e-7];

typerun    = 'u';  % user-provided inputs to SGP4 Matlab function
opsmode    = 'a';  % afspc approach ('air force space command')
whichconst =  72;  % WGS72 constants (radius, gravitational parameter)
arcsec2rad = pi / (180*3600);

spacecraft.ID = 28922;
spacecraft.Name = 'GIOVE-A';

satrec = read_TLE(spacecraft.ID,whichconst);

[year,mon,day,hr,minutes,sec] = invjday(satrec.jdsatepoch, satrec.jdsatepochf);
satrec_epoch_str = sprintf('%d-%02d-%02dT%02d:%02d:%02.6f', [year,mon,day,hr,minutes,sec]);
satrec_epoch_et = cspice_str2et(satrec_epoch_str);

% Nutation correction
ddpsi = -0.117144*arcsec2rad; %  [rad]
ddeps = -0.009050*arcsec2rad; %  [rad]
ateme = [0;0;0]; % teme acceleration

R = cspice_bodvrd('EARTH', 'RADII', 3);
flat = (R(1)-R(3))/R(1);

tformat = 'YYYY/MO/DD HR:MN:SC'; 

etf = cspice_str2et('2021-11-23T14:00:00');

N = floor((etf-satrec_epoch_et)/60);
et_vec = zeros(1,N+1);
et_vec(1) = satrec_epoch_et;
et = cspice_str2et('2021-11-19T14:26:00');
for i = 2:N+1
    et_vec(i) = et;
    et = et+60;
end

% Measurements recovery from: sgp4 propagation + random noise
for i = 1:length(station)
    if i == 1
        station(i).meas = zeros(length(station(i).visibility),3);
    else
        station(i).meas = zeros(length(station(i).visibility),2);
    end
    for j = 1:length(station(i).visibility)
        if station(i).visibility(j) == 1
            tsince = (et_vec(j) - satrec_epoch_et)/60.0; % minutes from TLE epoch
            [satrec,rteme_satrec,vteme_satrec] = sgp4(satrec,  tsince);
    
            % Compute centuries from 2000-01-01T00:00:00.00 TDT
            ttt = cspice_unitim(et_vec(i), 'ET', 'TDT')/cspice_jyear()/100;
            
            [reci_satrec, veci_satrec] = ...
                teme2eci(rteme_satrec, vteme_satrec, [0.0;0.0;0.0],  ttt, ddpsi, ddeps);
            if i == 1      % Milan, radar
                [Az, El, range, ~] = ...
                    antenna_pointing(station(i).Name, et_vec(j), [reci_satrec;veci_satrec]);
                station(i).meas(j,1:3) = normrnd([Az, El, range], [deg2rad(station(i).noise), deg2rad(station(i).noise), station(i).noise_range]);
            else     % Wellington / La-Silla, optical telescope
                reci_station = cspice_spkpos(station(i).Name, et_vec(j), 'J2000', 'NONE', 'EARTH');
                [~, Ra, Dec] = cspice_recrad(reci_satrec-reci_station);
                station(i).meas(j,1:2) = normrnd([Ra, Dec], [deg2rad(station(i).noise), deg2rad(station(i).noise)]);
            end
        end
    end
end

%%
% Batch filtering: LSM on variance without a priori info
etf = cspice_str2et('2021-11-21T14:00:00');
[~, idx] = min(abs(etf-et_vec));
t = et_vec(1:idx);
opt = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'iter');

% Keplerian motion
[Kep.x, Kep.resnorm, Kep.residual, Kep.exitflag, ~, ~, Kep.jac] = lsqnonlin(@(x) cost_fun(x, t, mu, station), initial.mean, [], [], opt);
Kep.jac = full(Kep.jac);
Kep.P_ls = Kep.resnorm/(length(Kep.residual)-length(initial.mean)).*inv(Kep.jac'*Kep.jac);

% J2 effect
j2 = 0.00108263;
[J2.x, J2.resnorm, J2.residual, J2.exitflag, ~, ~, J2.jac] = lsqnonlin(@(x) cost_fun(x, t, mu, station, j2), initial.mean, [], [], opt);
J2.jac = full(J2.jac);
J2.P_ls = J2.resnorm/(length(J2.residual)-length(initial.mean)).*inv(J2.jac'*J2.jac);


% % LSM on variance with a priori info
% % Keplerian motion
[Kep.x_api, Kep.resnorm_api, Kep.residual_api, Kep.exitflag_api, ~, ~, Kep.jac_api] = lsqnonlin(@(x) cost_fun_api(x, initial, t, mu, station), initial.mean, [], [], opt);
Kep.jac_api = full(Kep.jac_api);
Kep.P_ls_api = Kep.resnorm_api/(length(Kep.residual_api)-length(initial.mean)).*inv(Kep.jac_api'*Kep.jac_api);

% J2 effect
j2 = 0.00108263;
[J2.x_api, J2.resnorm_api, J2.residual_api, J2.exitflag_api, ~, ~, J2.jac_api] = lsqnonlin(@(x) cost_fun_api(x, initial, t, mu, station, j2), initial.mean, [], [], opt);
J2.jac_api = full(J2.jac_api);
J2.P_ls_api = J2.resnorm_api/(length(J2.residual_api)-length(initial.mean)).*inv(J2.jac_api'*J2.jac_api);


%% Additional propagation
% 2021-11-22T14:00:00
etf = cspice_str2et('2021-11-22T14:00:00');
[~, idx] = min(abs(etf-et_vec));
t = et_vec(1:idx);
opt = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'iter');

% Keplerian motion
[Kep.x, Kep.resnorm, Kep.residual, Kep.exitflag, ~, ~, Kep.jac] = lsqnonlin(@(x) cost_fun(x, t, mu, station), initial.mean, [], [], opt);
Kep.jac = full(Kep.jac);
Kep.P_ls = Kep.resnorm/(length(Kep.residual)-length(initial.mean)).*inv(Kep.jac'*Kep.jac);

% J2 effect
[J2.x, J2.resnorm, J2.residual, J2.exitflag, ~, ~, J2.jac] = lsqnonlin(@(x) cost_fun(x, t, mu, station, j2), initial.mean, [], [], opt);
J2.jac = full(J2.jac);
J2.P_ls = J2.resnorm/(length(J2.residual)-length(initial.mean)).*inv(J2.jac'*J2.jac);


% 2021-11-22T14:00:00
etf = cspice_str2et('2021-11-23T14:00:00');
[~, idx] = min(abs(etf-et_vec));
t = et_vec(1:idx);
opt = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'iter');

% Keplerian motion
[Kep.x, Kep.resnorm, Kep.residual, Kep.exitflag, ~, ~, Kep.jac] = lsqnonlin(@(x) cost_fun(x, t, mu, station), initial.mean, [], [], opt);
Kep.jac = full(Kep.jac);
Kep.P_ls = Kep.resnorm/(length(Kep.residual)-length(initial.mean)).*inv(Kep.jac'*Kep.jac);

% J2 effect
[J2.x, J2.resnorm, J2.residual, J2.exitflag, ~, ~, J2.jac] = lsqnonlin(@(x) cost_fun(x, t, mu, station, j2), initial.mean, [], [], opt);
J2.jac = full(J2.jac);
J2.P_ls = J2.resnorm/(length(J2.residual)-length(initial.mean)).*inv(J2.jac'*J2.jac);


%% Ex 3
clearvars; close all; clc; 

mu = cspice_bodvrd('SUN','GM', 3);

spacecraft.ID = '-143';
spacecraft.Name = 'exomars';

cspice_furnsh(['kernels\',spacecraft.Name, '.bsp']);

x_mean_0 = [  1.6847e8
             -1.0706e8
             -5.4724e7
              1.3435e1
              1.6871e1
              8.6608   ];

P_0 = [  2.01e4  -7.90    -4.05    -5.39e-3   6.37e-6  3.26e-6
        -7.90     2.01e4   2.64     6.37e-6  -5.38e-3 -2.07e-6
        -4.05     2.64     2.01e4   3.25e-6  -2.03e-6 -5.38e-3
        -5.39e-3  6.37e-6  3.25e-6  1.92e-7  -2.28e-9 -1.16e-9
         6.37e-6 -5.38e-3 -2.03e-6 -2.28e-9   1.91e-7  7.31e-10
         3.26e-6 -2.07e-6 -5.38e-3 -1.16e-9   7.31e-10 1.91e-7  ];

% Stations places
station(1).Name = 'Malargue';
station(2).Name = 'New-Norcia';

% Latitude
station(1).LAT = -35.77601;
station(2).LAT = -31.04823;

% Longitude
station(1).LON = -69.39819;
station(2).LON = 116.19147;

% Altitude
station(1).ALT = 1550; % m
station(2).ALT = 252; % m

% Type of sensor
station(1).type = 'Radar'; % monostatic
station(2).type = 'Radar';

% Noise
station(1).Az = deg2rad(1.5e-3);           station(1).El = deg2rad(1.3e-3);           station(1).range = 75e-3; % Az[°], El[°], range[km]
station(1).corr_Az_El = 0.1;               station(1).corr_Az_range = 0;              station(1).corr_El_range = 0;
% Noisy of 1st station
station(1).noise_mat = [                      station(1).Az^2                         station(1).corr_Az_El*station(1).Az*station(1).El        station(1).corr_Az_range*station(1).Az*station(1).range
                             station(1).corr_Az_El*station(1).Az*station(1).El                         station(1).El^2                         station(1).corr_El_range*station(1).El*station(1).range
                         station(1).corr_Az_range*station(1).Az*station(1).range    station(1).corr_El_range*station(1).El*station(1).range                        station(1).range^2                    ];

station(2).Az = deg2rad(1.5e-3);           station(2).El = deg2rad(1.3e-3);           station(2).range = 75e-3; % Az[°], El[°], range[km]
station(2).corr_Az_El = 0.1;               station(2).corr_Az_range = 0;              station(2).corr_El_range = 0;
% Noisy of 2nd station
station(2).noise_mat = [                      station(2).Az^2                         station(2).corr_Az_El*station(2).Az*station(1).El        station(2).corr_Az_range*station(2).Az*station(2).range
                             station(2).corr_Az_El*station(2).Az*station(2).El                         station(2).El^2                         station(2).corr_El_range*station(2).El*station(2).range
                         station(2).corr_Az_range*station(2).Az*station(2).range    station(2).corr_El_range*station(2).El*station(2).range                        station(2).range^2                    ];

% Minimum Elevation
station(1).min_El = 20; 
station(2).min_El = 20; 


% TDM data
directory = dir;
Files = dir(fullfile(directory(1).folder,'tdm'));
FileNames = cell(12,1);
for k = 3:length(Files)
    FileNames{k-2} = cellstr(Files(k).name);
    filePath = [ directory(1).folder,'\tdm\',char(FileNames{k-2}) ];
    tdm_table_new = read_TDM(filePath);
    if k == 3
        tdm_table = tdm_table_new;
    else
        tdm_table = [tdm_table; tdm_table_new];
    end
end
% sorting for reordering measures with chronological time
tdm_table = sortrows(tdm_table,2);


%%
% Unscented transform
alpha = 1.e-3;
beta = 2;
k = 0;
sz = size(tdm_table,1);

% Initialization 
x_hat_plus = zeros(sz-1,6);
P_plus = zeros(sz-1,6,6);

for j = 1:sz-1
    % Unscented weighted parameters
    lambda = alpha^2*(length(P_0)+k)-length(P_0);
    W0_m = lambda/(length(P_0)+lambda);
    Wi_m = 1/(2*(length(P_0)+lambda));
    W0_c = lambda/(length(P_0)+lambda)+(1-alpha^2+beta);
    Wi_c = Wi_m;
    sigmas = sigma_points_weighted(P_0, x_mean_0, lambda);
    sig_size = size(sigmas);
    station_name = tdm_table{j,1};
    meas = tdm_table{j,3:end};
    if strcmp(station_name,"MALARGUE") == 1
        R = inv(station(1).noise_mat);
    else
        R = inv(station(2).noise_mat);
    end
    % Unscented Filtering
    [y_hat, x_hat_minus, P_minus, P_ee, P_xy] = ukf(@(t,x) ode_orbit(t,x,mu), j, sigmas ,[W0_m; Wi_m*ones(sig_size(2)-1,1)], [W0_c; Wi_c*ones(sig_size(2)-1,1)], tdm_table, R); 

    % Update of mean and covariance
    K = P_xy*inv(P_ee);
    x_hat_plus(j,:) = x_hat_minus+K*(meas'-y_hat);
    P_plus(j,:,:) = P_minus-K*P_ee*K';

    % Restart cycle again
    x_mean_0 = x_hat_plus(j,:)';
    P_0 = squeeze(P_plus(j,:,:));
end

x_mean_final = x_mean_0;
P_final = P_0;


%% functions
% Ex3
function tdm = read_TDM(filename)

    fid = fopen(filename);

    % Read line 1-21
    for i = 1:21
        str = fgetl(fid);
        if i == 12
            new_str = split(str);
            station_name = new_str{3};
        end
    end

    str1 = fgetl(fid);
    count = 1;

    while strcmp(str1(1),'D') == 0

        % Read first line
        str2 = fgetl(fid);
        str3 = fgetl(fid);
        
        meas1 = split(str1);
        meas2 = split(str2);
        meas3 = split(str3);
        
        station(count,:) = station_name;
        Az(count) = str2double(meas1{4});
        t(count) = cspice_str2et(meas1{3});
        El(count) = str2double(meas2{4});
        Range(count) = str2double(meas3{4});

        str = fgetl(fid);
        str1 = fgetl(fid);
        count = count+1;

    end

    tdm = table(string(station), t', Az', El', Range');

    fclose(fid);

end

function [y_hat_minus, x_hat_minus, P_minus, P_ee, P_xy] = ukf(f, j, sigmas ,W_m, W_c, tdm_table, R)
%Unscented Transformation
%Input:
%        f: nonlinear map
%        t: time of the measurements
%        X: sigma points
%       Wm: weights for mean
%       Wc: weights for covraiance
%        n: numer of outputs of f
%        R: additive covariance
%Output:
%        y: transformed mean
%        Y: transformed sampling points
%        P: transformed covariance
%       Y1: transformed deviations

options_ode = odeset('reltol', 1e-8, 'abstol', 1e-8);
L = size(sigmas,2);

% Initialization
y_hat_minus = zeros(3,1);
dy = zeros(3,1);
x_hat_minus = zeros(6,1);
dx = zeros(6,1);
Y = zeros(3,L);
P_minus = zeros(6,6);
P_ee = zeros(3,3); 
P_xy = zeros(6,3);

station_name = char(tdm_table{j,1});
t(1) = tdm_table{j,2};
t(2) = tdm_table{j+1,2};
for k = 1:L
    % Sigma point propagation
    [~, x_ref] = ode113(f, t, sigmas(:,k), options_ode);
        
    reci_satrec = x_ref(end,1:3)';
    veci_satrec = x_ref(end,4:6)';

    % radar
    [Az, El, range] = ...
        antenna_pointing(station_name, t(2), [reci_satrec;veci_satrec]);

    % predicted measurements
    Y(:,k) = [Az; El; range];
    
    % UT weighted mean
    x_hat_minus = x_hat_minus+W_m(k)*sigmas(:,k);
    y_hat_minus = y_hat_minus+W_m(k)*Y(:,k);

end

for k = 1:L
    dx = dx+(sigmas(:,k)-x_hat_minus);
    dy = dy+(Y(:,k)-y_hat_minus); % array
end

for k = 1:L
    P_minus = P_minus+W_c(k).*dx*dx';
    P_ee = P_ee+W_c(k).*dy*dy'; 
    P_xy = P_xy+W_c(k).*dx*dy';
end

P_ee = P_ee+R;

end


% Ex2
function residual = cost_fun(x_mean_0, t, mu, station, varargin)
options_ode = odeset('reltol', 1e-10, 'abstol', 1e-10); % high tolerances just for better drawing
if length(varargin) == 1
    j2 = varargin{1};
    R =  cspice_bodvrd('EARTH', 'RADII', 3);
    [~, x_ref] = ode113(@(t,x) ode_orbit_J2(t,x,mu,j2,R(1)), t, x_mean_0, options_ode); 
else
    [~, x_ref] = ode113(@(t,x) ode_orbit(t,x,mu), t, x_mean_0, options_ode); 
end

% cost function
residual = [];

for i = 1:length(station)

    for j = 1:length(t)
        if station(i).visibility(j) == 1
            reci_satrec = x_ref(j,1:3)';
            veci_satrec = x_ref(j,4:6)';

            if i == 1      % Milan, radar
                [Az, El, range, ~] = ...
                    antenna_pointing(station(i).Name, t(j), [reci_satrec;veci_satrec]);

                % predicted measurements
                x_pred = [Az; El; range];
                W = diag([1/deg2rad(station(i).noise), 1/deg2rad(station(i).noise), 1/station(i).noise_range]);

            else     % Wellington / La-Silla, optical telescope
                reci_station = cspice_spkpos(station(i).Name, t(j), 'J2000', 'NONE', 'EARTH');
                [~, Ra, Dec] = cspice_recrad(reci_satrec-reci_station);

                % predicted measurements
                x_pred = [Ra; Dec];

                W = diag([1/deg2rad(station(i).noise), 1/deg2rad(station(i).noise)]);
            end

            % Update cost
            dy = x_pred-station(i).meas(j,:)';
            %         residual = residual+dy'*W*dy;
            residual = [residual; W*dy];

        end
    end
end

end

function residual = cost_fun_api(x_mean_0, initial, t, mu, station, varargin)
options_ode = odeset('reltol', 1e-10, 'abstol', 1e-10); % high tolerances just for better drawing
if length(varargin) == 1
    j2 = varargin{1};
    R =  cspice_bodvrd('EARTH', 'RADII', 3);
    [~, x_ref] = ode113(@(t,x) ode_orbit_J2(t,x,mu,j2,R(1)), t, x_mean_0, options_ode); 
else
    [~, x_ref] = ode113(@(t,x) ode_orbit(t,x,mu), t, x_mean_0, options_ode); 
end

% cost function
residual = [];

W_api = sqrtm(inv(initial.cov));
dy_api = W_api*(x_mean_0-initial.mean);
residual = [residual; dy_api];

for i = 1:length(station)
    for j = 1:length(t)
        if station(i).visibility(j) == 1   
            reci_satrec = x_ref(j,1:3)';
            veci_satrec = x_ref(j,4:6)';
    
            if i == 1      % Milan, radar
                [Az, El, range, ~] = ...
                    antenna_pointing(station(i).Name, t(j), [reci_satrec;veci_satrec]);
    
                % predicted measurements
                x_pred = [Az; El; range];
                W_m = diag([1/deg2rad(station(i).noise), 1/deg2rad(station(i).noise), 1/station(i).noise_range]);
    
            else     % Wellington / La-Silla, optical telescope
                reci_station = cspice_spkpos(station(i).Name, t(j), 'J2000', 'NONE', 'EARTH');
                [~, Ra, Dec] = cspice_recrad(reci_satrec-reci_station);
    
                % predicted measurements
                x_pred = [Ra; Dec];
    
                W_m = diag([1/deg2rad(station(i).noise), 1/deg2rad(station(i).noise)]);
            end
    
            % Update cost
            dy = x_pred-station(i).meas(j,:)';
    %         residual = residual+dy'*W*dy;
            residual = [residual; W_m*dy];
    
        end
    end
end

end

function dy = ode_orbit_J2( ~, y, mu, J2, Re)
%
%      ode_orbit_perturb.m - Two body problem equation
%     
%     PROTOTYPE:
%     	dy = ode_orbit_perturb( ~, y, mu)
%     
%     DESCRIPTION:
%       This function calculates the derivative of position and velocity of
%       the SC with time considering J2 perturbation. Use for integration.
%     
%     INPUT:
%       y[6,1]  Vector with position and velocity of the SC [km,km/s].
%       mu[1]   Planet constant [km3/s2]
%       J2[1]   Planet oblateness parameter.
%       Re[1]   Planet mean radius [km]
%     
%     OUTPUT:
%       dy[]    Vector with velocity and acceleration of the SC [km/s,km/s2].
%     
%     CALLED FUNCTIONS:
%       -
%     LAST UPDATED:
%      20/01/2020
%
%     CREATED BY:
%      Bardazzi N., Carcano S., Domaschio J., Maestro Martinez J.D.


    a_vector_x = (y(1)/norm(y(1:3)))*(5*y(3)^2/norm(y(1:3))^2-1);
    a_vector_y = (y(2)/norm(y(1:3)))*(5*y(3)^2/norm(y(1:3))^2-1);
    a_vector_z = (y(3)/norm(y(1:3)))*(5*y(3)^2/norm(y(1:3))^2-3);
    aJ2 = [a_vector_x;a_vector_y;a_vector_z];
    aJ2 = aJ2*((3/2)*(J2*mu*Re^2)/norm(y(1:3))^4);
    dy = [y(4:6);(-mu/norm(y(1:3))^3)*y(1:3)+aJ2];
end

% Ex1
function [X, Y, Z] = cylinder2P(R, N,r1,r2)
%  CYLINDER:  A function to draw a N-sided cylinder based on the
%             generator curve in the vector R.
%
%  Usage:      [X, Y, Z] = cylinder(R, N)
%
%  Arguments:  R - The vector of radii used to define the radius of
%                  the different segments of the cylinder.
%              N - The number of points around the circumference.
%
%  Returns:    X - The x-coordinates of each facet in the cylinder.
%              Y - The y-coordinates of each facet in the cylinder.
%              Z - The z-coordinates of each facet in the cylinder.
%
%  Author:     Luigi Barone
%  Date:       9 September 2001
%  Modified:   Per Sundqvist July 2004

    % The parametric surface will consist of a series of N-sided
    % polygons with successive radii given by the array R.
    % Z increases in equal sized steps from 0 to 1.

    % Set up an array of angles for the polygon.
    theta = linspace(0,2*pi,N);

    m = length(R);                 % Number of radius values
                                   % supplied.

    if m == 1                      % Only one radius value supplied.
        R = [R; R];                % Add a duplicate radius to make
        m = 2;                     % a cylinder.
    end


    X = zeros(m, N);             % Preallocate memory.
    Y = zeros(m, N);
    Z = zeros(m, N);
    
    v=(r2-r1)/sqrt((r2-r1)*(r2-r1)');    %Normalized vector;
    %cylinder axis described by: r(t)=r1+v*t for 0<t<1
    R2=rand(1,3);              %linear independent vector (of v)
    x2=v-R2/(R2*v');    %orthogonal vector to v
    x2=x2/sqrt(x2*x2');     %orthonormal vector to v
    x3=cross(v,x2);     %vector orthonormal to v and x2
    x3=x3/sqrt(x3*x3');
    
    r1x=r1(1);r1y=r1(2);r1z=r1(3);
    r2x=r2(1);r2y=r2(2);r2z=r2(3);
    vx=v(1);vy=v(2);vz=v(3);
    x2x=x2(1);x2y=x2(2);x2z=x2(3);
    x3x=x3(1);x3y=x3(2);x3z=x3(3);
    
    time=linspace(0,1,m);
    for j = 1 : m
      t=time(j);
      X(j, :) = r1x+(r2x-r1x)*t+R(j)*cos(theta)*x2x+R(j)*sin(theta)*x3x; 
      Y(j, :) = r1y+(r2y-r1y)*t+R(j)*cos(theta)*x2y+R(j)*sin(theta)*x3y; 
      Z(j, :) = r1z+(r2z-r1z)*t+R(j)*cos(theta)*x2z+R(j)*sin(theta)*x3z;
    end

    %surf(X, Y, Z);
end

function h = plot_gaussian_ellipsoid(m, C, sdwidth, npts, axh)
% PLOT_GAUSSIAN_ELLIPSOIDS plots 2-d and 3-d Gaussian distributions
%
% H = PLOT_GAUSSIAN_ELLIPSOIDS(M, C) plots the distribution specified by 
%  mean M and covariance C. The distribution is plotted as an ellipse (in 
%  2-d) or an ellipsoid (in 3-d).  By default, the distributions are 
%  plotted in the current axes. H is the graphics handle to the plotted 
%  ellipse or ellipsoid.
%
% PLOT_GAUSSIAN_ELLIPSOIDS(M, C, SD) uses SD as the standard deviation 
%  along the major and minor axes (larger SD => larger ellipse). By 
%  default, SD = 1. Note: 
%  * For 2-d distributions, SD=1.0 and SD=2.0 cover ~ 39% and 86% 
%     of the total probability mass, respectively. 
%  * For 3-d distributions, SD=1.0 and SD=2.0 cover ~ 19% and 73%
%     of the total probability mass, respectively.
%  
% PLOT_GAUSSIAN_ELLIPSOIDS(M, C, SD, NPTS) plots the ellipse or 
%  ellipsoid with a resolution of NPTS (ellipsoids are generated 
%  on an NPTS x NPTS mesh; see SPHERE for more details). By
%  default, NPTS = 50 for ellipses, and 20 for ellipsoids.
%
% PLOT_GAUSSIAN_ELLIPSOIDS(M, C, SD, NPTS, AX) adds the plot to the
%  axes specified by the axis handle AX.
%
% Examples: 
% -------------------------------------------
%  % Plot three 2-d Gaussians
%  figure; 
%  h1 = plot_gaussian_ellipsoid([1 1], [1 0.5; 0.5 1]);
%  h2 = plot_gaussian_ellipsoid([2 1.5], [1 -0.7; -0.7 1]);
%  h3 = plot_gaussian_ellipsoid([0 0], [1 0; 0 1]);
%  set(h2,'color','r'); 
%  set(h3,'color','g');
% 
%  % "Contour map" of a 2-d Gaussian
%  figure;
%  for sd = [0.3:0.4:4],
%    h = plot_gaussian_ellipsoid([0 0], [1 0.8; 0.8 1], sd);
%  end
%
%  % Plot three 3-d Gaussians
%  figure;
%  h1 = plot_gaussian_ellipsoid([1 1  0], [1 0.5 0.2; 0.5 1 0.4; 0.2 0.4 1]);
%  h2 = plot_gaussian_ellipsoid([1.5 1 .5], [1 -0.7 0.6; -0.7 1 0; 0.6 0 1]);
%  h3 = plot_gaussian_ellipsoid([1 2 2], [0.5 0 0; 0 0.5 0; 0 0 0.5]);
%  set(h2,'facealpha',0.6);
%  view(129,36); set(gca,'proj','perspective'); grid on; 
%  grid on; axis equal; axis tight;
% -------------------------------------------
% 
%  Gautam Vallabha, Sep-23-2007, Gautam.Vallabha@mathworks.com

%  Revision 1.0, Sep-23-2007
%    - File created
%  Revision 1.1, 26-Sep-2007
%    - NARGOUT==0 check added.
%    - Help added on NPTS for ellipsoids

if ~exist('sdwidth', 'var'), sdwidth = 1; end
if ~exist('npts', 'var'), npts = []; end
if ~exist('axh', 'var'), axh = gca; end

if numel(m) ~= length(m), 
    error('M must be a vector'); 
end
if ~( all(numel(m) == size(C)) )
    error('Dimensionality of M and C must match');
end
if ~(isscalar(axh) && ishandle(axh) && strcmp(get(axh,'type'), 'axes'))
    error('Invalid axes handle');
end

set(axh, 'nextplot', 'add');

switch numel(m)
   case 2, h=show2d(m(:),C,sdwidth,npts,axh);
   case 3, h=show3d(m(:),C,sdwidth,npts,axh);
   otherwise
      error('Unsupported dimensionality');
end

if nargout==0,
    clear h;
end

%-----------------------------
function h = show2d(means, C, sdwidth, npts, axh)
if isempty(npts), npts=50; end
% plot the gaussian fits
tt=linspace(0,2*pi,npts)';
x = cos(tt); y=sin(tt);
ap = [x(:) y(:)]';
[v,d]=eig(C); 
d = sdwidth * sqrt(d); % convert variance to sdwidth*sd
bp = (v*d*ap) + repmat(means, 1, size(ap,2)); 
h = plot(bp(1,:), bp(2,:), '-', 'parent', axh);
end

%-----------------------------
function h = show3d(means, C, sdwidth, npts, axh)
if isempty(npts), npts=20; end
[x,y,z] = sphere(npts);
ap = [x(:) y(:) z(:)]';
[v,d]=eig(C); 
if any(d(:) < 0)
   fprintf('warning: negative eigenvalues\n');
   d = max(d,0);
end
d = sdwidth * sqrt(d); % convert variance to sdwidth*sd
bp = (v*d*ap) + repmat(means, 1, size(ap,2)); 
xp = reshape(bp(1,:), size(x));
yp = reshape(bp(2,:), size(y));
zp = reshape(bp(3,:), size(z));
h = surf(axh, xp,yp,zp);
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
%      A_mat
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

function sigmas = sigma_points_weighted(P,mean,lambda)
sigmas = zeros(6,2*length(P)+1);
sqr_arg = sqrtm((length(P)+lambda)*P);
sigmas(:,1) = mean;
for i = 2:2:length(P)*2+1
    sigmas(:,i) = mean+sqr_arg(i/2,:)';
    sigmas(:,i+1) = mean-sqr_arg(i/2,:)';
end

end

function [sat_azimuth, sat_elevation, sat_range, sat_range_rate, ROT_ECI2TOPO] = antenna_pointing(stationName, et, rv_target_eci)

%station_visibility Summary of this function goes here
%   Detailed explanation goes here

topoFrame = [stationName, '_TOPO'];

% Get from SPK kernel the station state in ECI reference frame
rv_station_eci = cspice_spkezr(stationName, et, 'J2000', 'NONE', 'EARTH');

% Compute state vector of satellite as seen from the station in J2000
rv_station_sat_eci = rv_target_eci - rv_station_eci;

% Get state rotation matrix from ECI (J2000) to TOPOCENTRIC
ROT_ECI2TOPO = cspice_sxform('J2000',topoFrame, et);

% Convert target ECI state into topocentric
rv_station_sat_topo = zeros(size(rv_station_sat_eci));
for i = 1:size(rv_station_sat_eci,2)
    rv_station_sat_topo(:,i) = ROT_ECI2TOPO(:,:,i)*rv_station_sat_eci(:,i);
end

% Compute range and range-rate
% Compute euclidean norm along direction 1 (column-wise)
sat_range      = vecnorm(rv_station_sat_topo(1:3,:),2,1);

% Range rate is the scalar product of the velocity and a unit vector in the
% range direction
sat_range_rate = dot(rv_station_sat_topo(4:6,:), rv_station_sat_topo(1:3,:)./sat_range);

% Compute azimuth and elevation
rll_station_sat = cspice_xfmsta(rv_station_sat_topo,'RECTANGULAR','LATITUDINAL','EARTH');

%fprintf('%.5f %.5f\n',sat_range-rll_station_sat(1,:))
%fprintf('%.5f %.5f\n',sat_range_rate-rll_station_sat(4,:))

sat_azimuth = rll_station_sat(2,:);
sat_elevation = rll_station_sat(3,:);

end
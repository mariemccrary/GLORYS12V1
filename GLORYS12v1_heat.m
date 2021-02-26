% Marie McCrary Jan 2, 2021

% This file loads temperature data from monthly GLORYS12v1, subsets it
% by polygon mask, and calculates heat content within that mask.

% The first data file contans the data for each of the reanalyses. The
% second data file contains the mean and standard deviation of variables.

% Several calculations include functions from the Climate Data Toolbox:
% Chad A. Greene, Kaustubh Thirumalai, Kelly A. Kearney, José Miguel
% Delgado, Wolfgang Schwanghart, Natalie S. Wolfenbarger, Kristen M. Thyng,
% David E. Gwyther, Alex S. Gardner, and Donald D. Blankenship. The Climate
% Data Toolbox for MATLAB. Geochemistry, Geophysics, Geosystems 2019.
% doi:10.1029/2019GC008392 https://doi.org/10.1029/2019GC008392

% Bounding boxes were calculated by determining lat, lon coordinates from a
% polygon from https://boundingbox.klokantech.com/, and converting to i,j
% indeces 

% data = ('global-reanalysis-phy-001-031-grepv2-monthly_temp.nc');
data = ('CMEMS_global-reanalysis-phy-001-030-monthly_1993_2019.nc');

c = 4186; % J/kgC
rho = 1027; % kg/m3

% Full grid space
lats = double(ncread(data, 'latitude'));
lons = double(ncread(data, 'longitude'));
depth = double(ncread(data, 'depth'));

[lat, lon] = meshgrid(lats, lons);

% Temperature data (longitude,latitude,depth,time)

% t_oras = ncread(data, 'thetao_oras');
% t_cglo = ncread(data, 'thetao_cglo');
% t_glor = ncread(data, 'thetao_glor');
% t_foam = ncread(data, 'thetao_foam');

thetao = ncread(data, 'thetao', [1 1 1 1], [Inf Inf 1 Inf]); 
thetao = squeeze(thetao);

% Time variable
time = ncread(data, 'time'); % hours since 1950-01-01
timespan = length(time);

dtime = datetime(1950, 1, 1, time, 0, 0)

date = datevec(dtime);
time = datenum(date); 

%% Grid cell volume

A = cdtarea(lat, lon);

Vol = zeros(length(lons), length(lats), length(depth));
Vol_temp = zeros(length(lons), length(lats), length(depth));

for i = 1:length(depth)
    Vol(:,:,i) = A .* depth(i);
end

for i = 1:8
    Vol_temp(:,:,i) = Vol(:,:,i);
end
Vol10 = nansum(Vol_temp,3);
clear Vol_temp

%% Depth Avg Temperature

d = zeros(length(depth),1);

d(1) = 2*depth(1);
dtot = sum(d);

for i=2:25
     d(i) = 2*(depth(i)-dtot);
     dtot = sum(d);
end

Tsurf(:,:,1,:) = thetao(:,:,1,:);
Tsurf = squeeze(Tsurf);

% Depth Avg Temperature to 9.8 m, z=8

T10_mean = zeros(length(lons), length(lats), length(depth), timespan);
for i=1:8
    T10_mean(:,:,i,:) = (thetao(:,:,i,:).*d(i));
end
 
T10d = nansum(T10_mean,3); % Sums T along 3rd dim (depth)
T10d = squeeze(T10d);
Tbar10_mean = T10d./depth(8);

%% Heat content in each layer

% Heat is defined as Q = m*c*deltaT = rho*Vol*c*deltaT

Q_surf = zeros(length(lons), length(lats), length(time));
for i = 1:length(time)
    Q_surf(:,:,i) = rho.*c.* Vol(:,:,1) .*thetao(:,:,i);
end

Q_10 = zeros(length(lons), length(lats), timespan);

for i = 1:timespan
    Q_10(:,:,i) = rho.*c.* Vol10(:,:) .*Tbar10_mean(:,:,i);
end


%% Subset data by location

% Mackenzie Bounding box
bb_m = [-136.3858148812 69.0509501902;  -137.9436868227 69.3128192368 ;  -139.3367765793  69.6737303448 ;  -139.4190164954  70.4112187033 ;  -138.9257281747  71.0003029587 ;  -137.0775977055 71.4307629332 ;  -135.7256279071  71.3044102369 ;  -132.8894090333  71.0957989552 ;  -132.3757119098 70.4830637218 ;  -131.9398705894 69.8884664616 ;  -132.9863130522 69.6935711514 ;  -133.8059910402 69.5899640842 ;  -134.4672362486 69.7206254855 ;  -134.915474027 69.5695514246 ;  -135.304753042 69.4292489187 ;  -135.8775601052 69.2859396466 ;-136.3858148812 69.0509501902];


% 2007 trajectory

% Entire trajectory
bb = [-138.4553873539 69.6480360831; -168.739887476 80.0479196392; -164.3470728397 80.5206104623; -133.2729303837 70.5811562615; -138.4553873539 69.6480360831];
X = bb(:, 2);
Y = bb(:, 1);

% Section 1
bb1 = [-137.6319030466 69.4569141084; -148.3834397601 73.9334964563; -144.8224918643 74.7839890203; -133.713628199 70.6698161604; -137.6319030466 69.4569141084];
X1 = bb1(:, 2);
Y1 = bb1(:, 1);

% Section 2
bb2 = [-149.2669617035 73.5193651501; -161.7711197235 77.3049936832; -157.6148544647 78.0819441188; -145.3589187958 74.7150701086; -149.2669617035,73.5193651501];
X2 = bb2(:, 2);
Y2 = bb2(:, 1);

% Section 3
bb3 = [-161.7711197235 77.3049936832; -170.0579010346 79.9726420852; -165.859407077 80.5110516318; -157.6148544647 78.0819441188; -161.7711197235,77.3049936832]; 
X3 = bb3(:, 2);
Y3 = bb3(:, 1);

% Bering Strait bounding box
bb_b = [-168.2622826099 65.8531842878; -179.9361145496 69.5903310783; -163.8991653919 69.6070264563; -167.3826897144 68.5068058458; -163.8512718678 67.1302345154; -168.2622826099 65.8531842878];

% Kolyma bounding box
%bb_k = 

X_m = bb_m(:, 2);
Y_m = bb_m(:, 1);
boundm = geoshape(X_m, Y_m);
maskm = geomask(lat, lon, boundm.Latitude, boundm.Longitude);

X_b = bb_b(:, 2);
Y_b = bb_b(:, 1);
boundb = geoshape(X_b, Y_b);
maskb = geomask(lat, lon, boundb.Latitude, boundb.Longitude);

bound1 = geoshape(X1,Y1);
bound2 = geoshape(X2,Y2);
bound3 = geoshape(X3,Y3);
mask1 = geomask(lat, lon, bound1.Latitude, bound1.Longitude);
mask2 = geomask(lat, lon, bound2.Latitude, bound2.Longitude);
mask3 = geomask(lat, lon, bound3.Latitude, bound3.Longitude);

%% Plot bounding boxes
figure
hold on
worldmap([65 80], [-180 -130])
geoshow('landareas.shp', 'FaceColor', [0.85 0.85 0.85])
geoshow('worldrivers.shp','Color', 'blue')
geoshow(boundm.Latitude, boundm.Longitude, 'DisplayType', 'polygon', 'FaceColor', 'blue')
geoshow(boundb.Latitude, boundb.Longitude, 'DisplayType', 'polygon', 'FaceColor', 'red')
title('Bounding Boxes')

%% Heat Content of Mackenzie and Bering bounding boxes

% Mackenzie
Vol_surf =  Vol(:,:,1);
Vol_surf = squeeze(Vol_surf);

for i = 1:length(time)
    Q_surf(:,:,i) = rho.*c.* Vol_surf(:,:) .*thetao(:,:,i);
end

Tbar_surf_m = zeros(length(lons), length(lats), length(time));

for i=1:length(time)
    Tbar_surf_m(:,:,i) = thetao(:,:,i).*maskm(:,:);
end

Q_surf_m = zeros(length(lons), length(lats), length(time));


for i = 1:length(time)
    Q_surf_m(:,:,i) = rho.*c.* Vol_surf(:,:) .*Tbar_surf_m(:,:,i);
end

Q_surf_m_sum1 = nansum(Q_surf_m, [1 2]);
Q_surf_m_sum1 = squeeze(Q_surf_m_sum1);

% Bering Strait

Tbar10_b = zeros(length(lons), length(lats), timespan);

for i=1:timespan
    Tbar10_b(:,:,i) = Tbar10_mean(:,:,i).*maskb(:,:);
end

Q_10_b = zeros(length(lons), length(lats), timespan);


for i = 1:timespan
    Q_10_b(:,:,i) = rho.*c.* Vol10(:,:) .*Tbar10_b(:,:,i);
end

Q_10_b_sum1 = sum(Q_10_b, [1 2]);
Q_10_b_sum1 = squeeze(Q_10_b_sum1);

figure
hold on
plot(time, Q_surf_m_sum1, 'b')
%plot(time, Q_10_b_sum1, 'r')

datetick('x','mm/yy')
xlabel('time')
ylabel('Q(J)')
pm = polyfit(time, Q_surf_m_sum1, 1);
%pb = polyfit(time, Q_10_b_sum1, 1);
coefm = polyval(pm, time);
%coefb = polyval(pb, time);
plot(time, coefm, 'b.')
%plot(time, coefb, 'r.')

trend(Q_surf_m_sum1, time)

%%

Tbar10_M1 = zeros(length(lons), length(lats), timespan);
Tbar10_M2 = zeros(length(lons), length(lats), timespan);
Tbar10_M3 = zeros(length(lons), length(lats), timespan);

for i=1:timespan
    Tbar10_M1(:,:,i) = Tbar10_mean(:,:,i).*mask1(:,:);
end
for i=1:timespan
    Tbar10_M2(:,:,i) = Tbar10_mean(:,:,i).*mask2(:,:);
end
for i=1:timespan
    Tbar10_M3(:,:,i) = Tbar10_mean(:,:,i).*mask3(:,:);
end

Q_10_M1 = zeros(length(lons), length(lats), timespan);
Q_10_M2 = zeros(length(lons), length(lats), timespan);
Q_10_M3 = zeros(length(lons), length(lats), timespan);

for i = 1:timespan
    Q_10_M1(:,:,i) = rho.*c.* Vol10(:,:) .*Tbar10_M1(:,:,i);
end

Q_10_M_sum1 = sum(Q_10_M1, [1 2]);
Q_10_M_sum1 = squeeze(Q_10_M_sum1);

for i = 1:timespan
    Q_10_M2(:,:,i) = rho.*c.* Vol10(:,:) .*Tbar10_M2(:,:,i);
end

Q_10_M_sum2 = sum(Q_10_M2, [1 2]);
Q_10_M_sum2 = squeeze(Q_10_M_sum2);

for i = 1:timespan
    Q_10_M3(:,:,i) = rho.*c.* Vol10(:,:) .*Tbar10_M3(:,:,i);
end

Q_10_M_sum3 = sum(Q_10_M3, [1 2]);
Q_10_M_sum3 = squeeze(Q_10_M_sum3);

figure
hold on
plot(time, Q_10_M_sum1, 'b')
plot(time, Q_10_M_sum2, 'r')
plot(time, Q_10_M_sum3, 'k')
dateFormat = 11;
datetick('x','mm/yy')
xlabel('time')
ylabel('Q(J)')
p = polyfit(time, Q_10_M_sum, 1);
trend(Q_10_M_sum, time)
%% Plot full data and subsetted data to compare

figure
hold on
worldmap([65 80], [-180 -130])
pcolorm(lat, lon, (Q_10_M(:,:,57))./(10^17));
cmocean('thermal')
caxis([0 6])
h=colorbar;
h.Location = 'southoutside';
geoshow('landareas.shp', 'FaceColor', [0.85 0.85 0.85])
geoshow('worldrivers.shp','Color', 'blue')
ylabel(h,'J (x 10^{17})');
mlabel('off');
plabel('off');
title('10 m','FontWeight','normal')
%% NOTES %%
%Tbar, Htotal and 
%estimate, FWtotal, Sbar, Rhobar
% h=colorbar('SouthOutside');
% set(h, 'Position', [.1 .05 .8150 .05]);
% h.Label.String= 'Sea Ice Concentration (%)';
% %% Ice Data
 [ci,lat,lon] = arcticseaice('September 15, 2007', 'km', 'noplot');
% 
figure
hold on
worldmap([65 90], [-180 -130])
pcolorm(lat, lon, ci)
cmocean('ice')
h=colorbar;
caxis([0 100])
h.Location = 'southoutside';
geoshow('landareas.shp', 'FaceColor', [0.85 0.85 0.85])
geoshow('worldrivers.shp','Color', 'blue')
hold on
ylabel(h,'Concentration (%)')
mlabel('off');
plabel('off');
title('Ice','FontWeight','normal')
% 
% % savefig('')


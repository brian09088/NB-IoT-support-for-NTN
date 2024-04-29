clc,clear,close all

%% LEO satellite [starlink-1007]
% 衛星目錄訊號(NORAD-ID) : 44713
% NSSDC ID : 2019-074A

%% Time Interval / Create and Visualize Scenario
startTime = datetime(2024,3,15,12,0,0);
stopTime = startTime + minutes(30);
sampleTime = 30; % seconds
sc = satelliteScenario(startTime,stopTime,sampleTime);
viewer = satelliteScenarioViewer(sc,ShowDetails=false);

%% user defined data : Orbit
% orbit propagator
OrbitPropagator = "sgp4"; % two−body−keplerian, sgp4 or sdp4
GSpointing = "satellite"; % zenith ( the antenna of the gs points to zenith ) or satellite (points to the satellite)        
Phasingmethod = "alternate"; % standard or alternative

%% Add Satellite Constellation
% As of February 2024, there are 5,438 Starlink satellites in orbit
numSat_Per_Plane = 11;
numOrbits = 6;
% defined satellite constellation
% ex: space-X starlink satellite dconstellation(satellite_tle.csv)
sat_tle = uigetfile('*.csv','Choose a .csv file');
starlink = satellite(sc, sat_tle,"OrbitPropagator", OrbitPropagator);

%% user defined data : Satellite
Re = 6371e3;          % Radius of Earth(km)
h0 = 554e3;           % altitude of satellite ex:starlink-1007
rP = Re + h0;         % Perigee of orbit (m) 近地點軌道高度
rA = Re + h0;         % Apogee of orbit (m) 遠地點軌道高度
semi_major_axis = (rA + rP)/2;  % 衛星軌道半長軸

%% user defind data : Satellite location
% get from defined satellite constellation (sat_tle)
% starlink_tle_1007 Launch date: November 11, 2019

%% user defind data : Ground Station locaiton
gs_name = '國家太空中心-TASA';
gs_lat = 24.4703;
gs_lon = 121.0015;
gs_location = [gs_lat, gs_lon];

% 圓心角(天北極夾角) contact angle (degree)
phi = calculat_phi(gs_lat, gs_lon);    

alpha = 53.1;             % inclination 發射傾角 elevation angle (>= 10 degree)
eccentricity = 0;         % 偏心率 e = 0.5*(rA - rP)/a
arg_of_periapsis = 0;     % 近日點幅角 degrees
trueAnomaly = 0;          % 真近點角(衛星與近地點之間的橢球焦點角距) degrees
raan = 0;                 % 升交點黃經 Right ascension of ascending node, in degrees
% 角度轉弧度 rad = deg * deg2rad

%% LEO satellites parameters

% transmitting antenna power 23(dBm) = 200(mW)
antenna_power_dBm = 23;     % P = 10^(dBm/10) -> 10*log(200) = 23
antenna_power_mW = dBm_to_mW(antenna_power_dBm);

sat_antenna_gain = 30;  % Antenna gain = 0 (dBi)

sat_aperture = 2;       % 天線孔徑(m)
antenna_gain_to_noise = 1.1;  % Antenna gain-to-noise-temperature (dB/K)

% antenna gain 天線增益
% When non-isotropic antennas are used(當採用非均向天線 ex: omnidirectional)
Gt = 0;                 % transmit antenna gain (dBi)
Gr = 0;                 % receive antenna gain (dBi)

EIRP = calculate_EIRP(antenna_power_dBm,sat_antenna_gain);

NF = 7;                 % noise figure (dB)
To = 290;               % 室溫/環境溫度 ambient temperature (K)
Ta = 290;               % antenna temperature (K)

%% Calculte free-space path loss
% frequency band of 2100 Hz for the center frequency of 2 GHz
signal_freq = 2.0e9; % signal frequency band : S-band (Hz) 
d = calculate_d(Re, alpha, h0);   % distance between sat & UE (meters)
FSPL = FSPL_model(signal_freq, d);

%% doppler shift pre-compensation
% Maximum Doppler shift 24 [ppm]

% Residual Doppler shift after pre-compensation
% 1.05 [ppm] for 50 km beam diameter
% 1.88 [ppm] for 90 km beam diameter
% 15.82 [ppm] for 1000 km beam diameter


LEO_Sat_1 = satellite(sc, semi_major_axis, ...
    eccentricity, ...
    inclination, ...
    raan, ...
    arg_of_periapsis, ...
    trueAnomaly, ...
    Name = "LEO 1", ...
    OrbitPropagator = "two-body-keplerian");

LEO_Sat_2 = satellite(sc,semi_major_axis, ...
    eccentricity, ...
    inclination, ...
    raan, ...
    arg_of_periapsis, ...
    trueanomaly, ...
    Name="LEO 2", ...
    OrbitPropagator="two-body-keplerian");

%% Add Grid of Ground Stations Covering Taiwan-ROC(地面基地台)
latlim = [22.00417 25.12825];
lonlim = [118.31833 121.753];
spacingInLatLon = 1; % 基站間距離，越小覆蓋率越高但計算時間較長(單位degrees)
EPSG_code = 3826; % (TWD97-EPSG:3826/3830)
proj = projcrs(EPSG_code); % projected coordinate reference system
spacingInXY = deg2km(spacingInLatLon)*1000; % 相對於spacingInLatLon轉換為xy座標系統(單位:meters)
[xlim,ylim] = projfwd(proj,latlim,lonlim);
R = maprefpostings(xlim,ylim,spacingInXY,spacingInXY);
[X,Y] = worldGrid(R);
[gridlat,gridlon] = projinv(proj,X,Y);

TW_landareas = readgeotable("E:\MATLAB\碩士論文\Coverage_Maps_for_Satellite_Constellation\gadm36_TWN_shapread\gadm36_TWN_0.shp");
% ROC = landareas(landareas.Name == "Taiwan",:);

T = geotable2table(TW_landareas,["Latitude","Longitude"]);
[landlat,landlon] = polyjoin(T.Latitude,T.Longitude);

bufwidth = 1;
[landlatb,landlonb] = bufferm(landlat,landlon,bufwidth,"outPlusInterior");
TW_liab = geopolyshape(landlatb,landlonb);

gridpts = geopointshape(gridlat,gridlon);
inregion = isinterior(TW_liab,gridpts);
gslat = gridlat(inregion);
gs_lon = gridlon(inregion);

gs = groundStation(sc,gslat,gs_lon);

%% antenna type 
% satellite : omnidirectional 全向天線
% ground stations : phased array antenna 相位陣列天線(iridium : custom 48-beam)

fq = 1e9; % Hz
txpower = 30; % dBW
beamWidth = 5; % degrees

% antennatype = custom 48-beam
if antennaType == "phasedarray" 
    tx = transmitter(sats, ...
        Frequency=fq, ...
        Power=txpower); 
    isotropic = arrayConfig(Size=[1 1]);
    isotropicAntenna(tx,Array= isotropic);
end

% Custom-48 Beam 天線：
if antennaType == "Custom 48-Beam"
    antenna = helperCustom48BeamAntenna(fq);
    tx = transmitter(sats, ...
        Frequency=fq, ...
        MountingAngles=[0,-90,0], ... % [yaw, pitch, roll] with -90 using Phased Array System Toolbox convention
        Power=txpower, ...
        Antenna=antenna);  
end

% Add Ground Station Receivers (地面接收站獲取訊號後透過基地台轉傳)
isotropic = arrayConfig(Size=[1 1]);
rx = receiver(gs,Antenna=isotropic);
pattern(tx,Size=500000);

% Compute Raster Coverage Map Data
delete(viewer)
maxsigstrength = satcoverage(gridpts,sc,startTime,inregion,beamWidth);

% 設定當地重要城市經緯度，進行定位與標註
TP_location = [25.02 121.31]; % Taipei
TC_location = [24.08 120.41]; % Taichung
Kao_location = [22.36 120.18]; % Kaosiung

textm(TP_location(1),TP_location(2),"Taipei")
textm(TC_location(1),TC_location(2),"Taichung")
textm(Kao_location(1),Kao_location(2),"Kaohsiung")


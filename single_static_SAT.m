% function for simulate satellite scenario

%% Static Satellite : (Geosynchronous)
% satellite : 1 (AsiaSat_9)
% GS : 1

Ground_Station();

% Create and Visualize Scenario
startTime = datetime(2024,4,11,12,0,0);
stopTime = startTime + hours(1);
sampleTime = 60; % seconds
sc = satelliteScenario(startTime,stopTime,sampleTime);
viewer = satelliteScenarioViewer(sc,ShowDetails=false);

OrbitPropagator = "two-body-keplerian"; % two-body-keplerian, sgp4 or sdp4
GSpointing = "satellite"; % zenith ( the antenna of the gs points to zenith ) or satellite (points to the satellite)

%% Use GEO-SAT : AsiaSAT_9 (Communications : TV)
% Add Satellites to the Satellite Scenario
tleFile = "AsiaSAT_9.csv";
sat = satellite(sc, tleFile, "OrbitPropagator", OrbitPropagator, Name="GEO_SAT");

% 設定當地重要城市經緯度，進行定位與標註
TP_location = [25.02 121.31]; % Taipei
TC_location = [24.08 120.41]; % Taichung
Kao_location = [22.36 120.18]; % Kaosiung

% Add Grid of Ground Stations Covering Taiwan-ROC(地面基地台)
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
gslon = gridlon(inregion);

% Add a camera to the GEO satellite
cam = conicalSensor(sat, 'Name', 'GEO Satellite Camera', 'MaxViewAngle', 90);

% Define the Geographical Site for the Ground Station
minElevationAngle = 30; % degrees
geoSite = groundStation(sc, ...
    'Name', gs_name, ...
    'Latitude', gs_lat, ...
    'Longitude', gs_lon, ...
    'MinElevationAngle', minElevationAngle);

% Visualize the Scenario
v = satelliteScenarioViewer(sc, 'ShowDetails', false);
sat(1).ShowLabel = true;
geoSite.ShowLabel = true;

% Visualize the Field Of View of the Camera
fov = fieldOfView(cam);
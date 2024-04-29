% function for simulate satellite scenario

%% Dynamic Satellite : 
% satellite : multi / LEO-constellation
% GS : 1

Ground_Station();

%% main function
startTime = datetime(2024,4,10,13,0,0);
stopTime = startTime + hours(6);
sampleTime = 30; % seconds
sc = satelliteScenario(startTime,stopTime,sampleTime);

OrbitPropagator = "sgp4"; % two−body−keplerian, sgp4 or sdp4
GSpointing = "satellite"; % zenith ( the antenna of the gs points to zenith ) or satellite (points to the satellite)        
Phasingmethod = "alternate"; % standard or alternative

% Add Satellites to the Satellite Scenario
tleFile = "leoSatelliteConstellation.tle";
sat = satellite(sc,tleFile, "OrbitPropagator", OrbitPropagator);

% Add Cameras to the Satellites
names = sat.Name + " Camera";
cam = conicalSensor(sat,"Name",names,"MaxViewAngle",90);

% Define the Geographical Site to be Photographed in the Satellite Scenario
minElevationAngle = 30; % degrees
geoSite = groundStation(sc, ...
    "Name",gs_name, ...
    "Latitude",gs_lat, ...
    "Longitude",gs_lon, ...
    "MinElevationAngle",minElevationAngle);

% Add Access Analysis Between the Cameras and the Geographical Site
ac = access(cam,geoSite);
ac.LineColor = 'red';

% Properties of access analysis objects
ac(1)

% Visualize the Scenario
v = satelliteScenarioViewer(sc,"ShowDetails",false);
sat(14).ShowLabel = true;
geoSite.ShowLabel = true;
show(sat(14));

% Visualize the Field Of View of the Camera
fov = fieldOfView(cam([cam.Name] == "Satellite 1 Camera"));
% Create and Visualize Scenario
startTime = datetime(2024,4,11,12,0,0);
stopTime = startTime + hours(1);
sampleTime = 60; % seconds
sc = satelliteScenario(startTime,stopTime,sampleTime);
viewer = satelliteScenarioViewer(sc, ShowDetails=false);

% Define the boundaries of Taiwan
taiwanBoundaries = [...
    20.5170, 118.2008; % Southwest corner
    25.3000, 122.0000; % Northwest corner
    25.3000, 121.3000; % Northeast corner
    20.6834, 119.6000; % Southeast corner
    20.5170, 118.2008  % Repeat the first point to close the polygon
]; 

% Create a geoshape object from the boundary coordinates
taiwanShape = geoshape(taiwanBoundaries(:, 1), taiwanBoundaries(:, 2));

% Calculate the centroid of the Taiwan boundaries
centroid = [mean(taiwanShape.Latitude), mean(taiwanShape.Longitude)];

% Calculate the azimuth angle from the satellite to the centroid of Taiwan
% Assuming the satellite is directly overhead
satelliteAzimuth = 180;

% Calculate the true anomaly based on the azimuth angle
trueanomaly = 90 - satelliteAzimuth;  % Assuming 90 degrees is directly overhead

% Update the geostationary satellite parameters with the adjusted true anomaly
semimajoraxis = 35784e3; % Semi-major axis for geostationary orbit (approximately)
inclination = 0; % Inclination for geostationary orbit
eccentricity = 0; % Eccentricity for geostationary orbit
argofperiapsis = 0; % Argument of periapsis for geostationary orbit
RAAN = 0; % Right ascension of ascending node for geostationary orbit

% Create the geostationary satellite with the adjusted true anomaly
satelliteName = 'Geostationary Satellite';
sats = satellite(sc, semimajoraxis, eccentricity, inclination, RAAN, argofperiapsis, trueanomaly, 'Name', satelliteName);

% Visualize the satellite scenario
v = satelliteScenarioViewer(sc, 'ShowDetails', false);
sat(1).ShowLabel = true;

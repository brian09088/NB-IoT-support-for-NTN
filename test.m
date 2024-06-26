clc,clear,close all

%%%%%% 地球同步衛星 %%%%%%

%% Time Interval / Create and Vddisualize Scenario
startTime = datetime(2024,4,2,5,0,0);
stopTime = startTime + hours(9);
sampleTime = 30; % seconds
sc = satelliteScenario(startTime,stopTime,sampleTime);
viewer = satelliteScenarioViewer(sc,ShowDetails=false);

%% User defined data - orbit propagator
% orbit propagator
OrbitPropagator = "sgp4"; % two−body−keplerian, sgp4 or sdp4
GSpointing = "satellite"; % zenith ( the antenna of the gs points to zenith ) or satellite (points to the satellite)        
Phasingmethod = "alternate"; % standard or alternative

%% Add Satellite Constellation

%% Add Satellite Constellation

% As of February 2024, there are 5,438 Starlink satellites in orbit
numSat_Per_Plane = 11;
numOrbits = 6;

% defined satellite constellation
% ex: space-X starlink satellite constellation(satellite_tle.csv)
sat_tle = uigetfile('*.csv','Choose a .csv file');

% Read the TLE file
data = readtable(sat_tle, 'ReadVariableNames', false, 'Delimiter', '\n');

% Determine the number of satellites in the TLE file
num_satellites = height(data) / 3;

% Create the satellite scenario
starlink_scenario = satelliteScenario;

% Add satellites to the scenario
for i = 1:num_satellites
    sat_name = data{3 * (i - 1) + 1, 1}{1};
    sat_line1 = data{3 * (i - 1) + 2, 1}{1};
    sat_line2 = data{3 * (i - 1) + 3, 1}{1};
    
    starlink_scenario.addSatellite(sat_name, sat_line1, sat_line2, "OrbitPropagator", OrbitPropagator);
end

% Show the satellite scenario
show(starlink_scenario);

%% User defined data - satellite parameters
Re = 6371e3;                    % Radius of Earth
h0 = 650e3;                     % altitude of satellite
LEO_altitude = h0;
rP = Re + LEO_altitude;         % Perigee of orbit (m) 近地點軌道高度
rA = Re + LEO_altitude;         % Apogee of orbit (m) 遠地點軌道高度
semi_major_axis = (rA + rP)/2;  % 衛星軌道半長軸
inclination = 53.1;            % 發射傾角 degrees
eccentricity = 0;               % 偏心率 e = 0.5*(rA - rP)/a
arg_of_periapsis = 0;           % 近日點幅角 degrees
trueAnomaly = 0;                % 真近點角(衛星與近地點之間的橢球焦點角距) degrees
raan = 0;                  % 升交點黃經 Right ascension of ascending node, in degrees
% 角度轉弧度 rad = deg * deg2rad

%% User defined data - satellites velocity (receiver speed)
v = 0;   % geosynchronous satellite : static (relative velovity = 0)
% v = ? ; % sun-geosynchronous satellite : dynamic (relative velovity = ?)

%% User defined data - UE parameters (support NB-IoT)
% Frequency band (fc) S-band (2 GHz)
% Antenna type and configuration (1,1,2) with omnidirectional antenna element
% Polarization Linear: +/− 45◦X − pol
% Antenna temperature = 290 [K]
% Noise figure = 7 [dB]
% Transmit power = 200 [mW]
% Antenna gain = 0 [dBi]

%% NB-IoT NTN Assumptions
% (i) The target UE is located within the spot beam
% (ii) the satellite can steer beams towards fixed points on earth using beamforming techniques
% (iii) assuming that the feeder and the inter-satellite link are ideal, the service link performance is analyzed
% (vi) a minimum elevation angle of 10 degrees is considered for the UE and the satellite

%% User defined data - beam frequency & coverage


%% User defined data - antenna type and beam coverage
% satellite : omnidirectional 全向天線(不須考慮衛星姿態控制，只需要設定軌道與飛行平面)
% ground stations : phased array antenna 相位陣列天線(iridium : custom 48-beam)

% **Regulatory Considerations** 
% S-band (2 to 4 GHz) 
% commonly used for satellite communication has allocated frequencies for such purposes by regulatory bodies like ITU
% Researchers may choose this band because of its availability and regulatory approval for satellite communication systems

% fq represents the frequency of the communication signal in Hertz (Hz)
fc = 2e9;       % center frequency in Hz(2GHz -> Hz)
fc_GHz = fc / (10^9);
fprintf("center frequency = %d GHz\n", fc_GHz)

% txpower represents the transmit power of the transmitter in dBW (decibels relative to 1 watt)
txpower_dBW = 30;   % dBW (dBm=dBW+30)
txpower_mW = dBW_to_mW(txpower_dBW);
fprintf("transmit power = %d dBW\n", txpower_dBW)

% beam_width represents the beam width of the phased array antenna in degrees
% This parameter determines the angular width of the main lobe of each beam
% 相控陣天線的波束寬度（degree）決定每個波束主瓣的角寬度
beam_width = 5; % degrees

antennaType = "Custom 48-Beam";
% if antennaType == "Custom 48-Beam" 
%   tx = transmitter(sats, ...
%        Frequency=fq, ...
%        Power=txpower); 
%        isotropic = arrayConfig(Size=[1 1]);
%        isotropicAntenna(tx,Array= isotropic);
% end

% UE_antenna_gain 裝置天線增益
UE_antenna_gain = 0;    % Antenna gain = 0 [dBi]
% Sat_antenna_gain 衛星天線增益
Sat_antenna_gain = 30;

% Transmit antenna gain 傳輸過程天線增益
% Transmit_antenna_gain = UE_antenna_gain + Sat_antenna_gain;
Transmit_antenna_gain = 0;

fprintf("UE_antenna_gain = %d \n", UE_antenna_gain)
fprintf("Sat_antenna_gain = %d \n", Sat_antenna_gain)
fprintf("Transmit_antenna_gain = %d \n", Transmit_antenna_gain)

% Transmit power = 200 [mW] (paper parameters setup)
UE_antenna_power = 200;  %(mW)
UE_antenna_dBm = mW_to_dBm(UE_antenna_power);
fprintf("UE_antenna_power = %f mW\n", UE_antenna_power)
fprintf("UE_antenna_power = %f dBm\n", UE_antenna_dBm)

% Satellite effective isotropic radiated power (EIRP) [dBW/MHz]
EIRP = calculate_EIRP(UE_antenna_power, UE_antenna_gain, Sat_antenna_gain);
fprintf("有效全向輻射功率 EIRP = %f (dBW)\n", EIRP) % -6.99

%% User defined data - Ground stations locaiton

gs_name = '國家太空中心-TASA';
gs_lat = 24.4703;
gs_lon = 121.0015;
gs_location = [gs_lat, gs_lon];

alpha = inclination;       % inclination (elevation angle)發射傾角 degrees (>10 degree)
% phi = calculate_phi(gs_lat, gs_lon); % central angle (地心_天北極 & 地心與衛星連線夾角)
% distance between sat & UE (meters)
d = calculate_d(Re, alpha, h0);
fprintf("communication path distance = %f (m)\n", d)

%% Maximum Doppler shift and the residual Doppler shift after pre-compensation
% Frequency band S-band (2 GHz)
% Satellite altitude 600 [km]
% Maximum Doppler shift 24 [ppm]

% Residual Doppler shift after precompensation
% 1.05 [ppm] for 50 km beam diameter
% 1.88 [ppm] for 90 km beam diameter
% 15.82 [ppm] for 1000 km beam diameter

%% Link-level simulation parameters
% Transmit power −6.99 dBW 200 mW ->10log(0.2)
% Transmit antenna gain 0 dBi Antenna gain
% EIRP (effective isotropic radiated power) −6.99 dBW Transmit power + antenna gain
% Atmospheric loss 0.07 dB [21]
% Shadow fading margin 3.00 dB
% Scintillation loss 2.20 dB
% Polarization loss 0.00 dB
% Additional losses 0.00 dB
% Boltzmann%s constant [k] −228.6 dBW/K/Hz
% Bandwidth [B] 52.55 dBHz

%% User defined data - Link Budget

% This is the satellite antenna
Tant.Up = 290; % Antenna noise temperature in Kelvin for uplink
% This is the ground station antenna
Tant.Down = 290; % Antenna noise temperature in Kelvin for downlink
Tant.Inter = 10; % Antenna noise temperature in Kelvin for interlink

NF = 2; % Noise figure in dB
T0 = 280; % Reference temperature in K
Pmin = -134; % Receiver sensitivity in dBW
freq.Inter = 2.2e9; % Frequency of interlink in Hz
freq.UpDown = 2.2e9; % Frequency of uplink/downlink in Hz
P_TR.Sat = 1; % Power in watts of the satellites
P_TR.GS = 5; % Power in watts of the ground stations
Bandwidth = 750e3; % Bandwidth in Hz
BitRate = 512e3; % Bit Rate in Hz (bit/s)

%% Calculte free-space path loss
FSPL = FSPL_model(fc, d);
fprintf("自由空間路徑損耗 FSPL = %f (dB)\n", FSPL)

%% Doppler shift model
fq_ref = fc;  % f reference
fq_edge = fc; % f edge
fq_res = fc;  % f residual
fd = Doppler_model(fc, v, alpha);
fprintf("Doppler model : %f (kHz)\n", fd)

%% Uplink Modulation Format & Order

% modulation format : SC-FDMA
% modulation order : QPSK/OQPSK
mod_order = 'OQPSK';

% SNR method
SNR_method = "berspec";

% BER specification
mod = mod_order;
% Bit Error Rate specification
ber_spec = 1e-6;

% Offset-Quadrature-Phase-Shift-Keying (OQPSK) is a variant of the QPSK modulation scheme 
% where the phase or timing of either the in-phase or 
% Quadrature component is shifted relative to each other by a one bit-period or 
% half a symbol-period Ts
% Modulation(mod_order , ber_spec)

EbNo_spec = Modulation(mod, ber_spec) ;
SNR_spec = EbNo_To_SNR(EbNo_spec, BitRate, Bandwidth) ;

%% Fading channel model 
% 1. ✅ 3GPP TDL-D
% 2. AWGN
fprintf("Fading channel model : 3GPP TDL-D \n")

% use MATLAB NTN_Channel_Model, call function file NTN_Channel_Model.m
% include HelperGenerateNTNChannel.m
% include HelperSetupNTNChannel.m
% include hNPDSCHInfo.m
NTN_Channel_Model();

%% Ch3 Problem Statement (Doppler compensation)

%% 3.1 Residual doppler shift : Doppler shift compensation & Fast Fourier Transform
FFT_size = 128;    % Fast Fourier Transform size (n=128)
n = 10;            % number of symbols between the DMRS symbols
Lcp = num2cell([10, repmat(9, 1, n-1)]);  % Assign 10 to the first symbol and 9 to the rest

% Number of time samples between channel estimation symbols
% L = calculate_L(Lcp, n, FFT_size);   
L = 960;

% Phase difference between channel estimation symbols, theta = 0
theta = calculate_theta();  

dopp_range = dopp_compensation(FFT_size, theta, L);  % range of the Doppler shift compensation

% maximum phase difference between DMRS symbols = pi, dopp_range = 1/15
% SCS = 15 kHz
SCS = 15;
fd_max = SCS * dopp_range;

% fd_max = 1kHz
fprintf("maximum possible Doppler shift compensation (fd_max) = %d (kHz)\n", fd_max);

%% Link-level simulation parameters
% FFT_size = 128;
SCS = 15;     %(kHz)
% Multiple access = SC-FDMA
% Modulation QPSK
% Channel coding 1/3 Turbo code
TBS = 120;    %(bits)
% Fading channel model 3GPP TDL-D
% Residual Doppler shift 0, 950 [Hz]
NRep = 2;
% Bandwidth = 180;    %(B = 180kHz)

%% 3.2 Link performance
% HARQ Configuration: 
% Configure parameters related to Hybrid Automatic Repeat reQuest (HARQ) 
% such as the number of retransmissions, coding scheme, redundancy versions 
% HARQ plays a crucial role in improving link reliability and performance

% Define Communication Scenarios: Define various communication scenarios considering different numbers of LEO satellites
% communication time variations, RTT (Round-Trip Time), packet sizes, and data transmission rates
% This will help in assessing the impact of satellite communication on link performance

%% Ch4 Proposed solution
% 4.1 Reduction in  beam coverage
% 4.2 Addition of DMRS symbol

%% Relay Algorithm 接下來進度
% Relay_algorithm()

%% 自定義衛星
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
    trueAnomaly, ...
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
gslon = gridlon(inregion);

gs = groundStation(sc,gslat,gslon);

% Add Ground Station Receivers (地面接收站獲取訊號後透過基地台轉傳)
isotropic = arrayConfig(Size=[1 1]);
tx = receiver(gs,Antenna=isotropic);


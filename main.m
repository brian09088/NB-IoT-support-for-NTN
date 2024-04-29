clc,clear,close all

%%%%%% 地球同步衛星 %%%%%%

%% Time Interval / Create and Vddisualize Scenario
% IRIS-A (NORAD-ID:51044) 
% Launch date: January 13, 2022
startTime = datetime(2024,4,20,8,0,0);
stopTime = startTime + hours(8.5);
sampleTime = 30; % seconds
sc = satelliteScenario(startTime,stopTime,sampleTime);
viewer = satelliteScenarioViewer(sc,ShowDetails=false);

%% User defined data - orbit propagator
% orbit propagator
OrbitPropagator = "sgp4"; % two−body−keplerian, sgp4 or sdp4
GSpointing = "satellite"; % zenith ( the antenna of the gs points to zenith ) or satellite (points to the satellite)        
Phasingmethod = "alternate"; % standard or alternative

%% Add Satellite Constellation
% As of February 2024, there are 5,438 Starlink satellites in orbit
numSat_Per_Plane = 11;
numOrbits = 6;
% defined satellite constellation
% ex: space-X Starlink satellite constellation(satellite_tle.csv)
sat_tle = uigetfile('*.csv','Choose a .csv file');
% sat_tle = tleread('***.tle');
starlink = satellite(sc, sat_tle, "OrbitPropagator",OrbitPropagator);

%% User defined data - satellite parameters
Re = 6371e3;                    % Radius of Earth (m)
h0 = 650e3;                     % altitude of satellite (m)
% IRIS-A : 486 km
% paper starlink : 552 km
SAT_altitude = h0;
rP = Re + SAT_altitude;         % Perigee of orbit (m) 近地點軌道高度
rA = Re + SAT_altitude;         % Apogee of orbit (m) 遠地點軌道高度
semi_major_axis = (rA + rP)/2;  % 衛星軌道半長軸
inclination = 53.1;             % 發射傾角 degrees (elevation angle) 53.1
eccentricity = 0;               % 偏心率 e = 0.5*(rA - rP)/a
arg_of_periapsis = 0;           % 近日點幅角 degrees
trueAnomaly = 0;                % 真近點角(衛星與近地點之間的橢球焦點角距) degrees
raan = 0;                  % 升交點黃經 Right ascension of ascending node, in degrees
% 角度轉弧度 rad = deg * deg2rad
h0_km = h0 / 1e3;
fprintf("Satellite altitude = %d km\n",h0_km)

%% User defined data - satellites velocity (receiver speed)
v = 0;   % geosynchronous satellite : static (relative velovity = 0)
% v = ? ; % sun-geosynchronous satellite : dynamic (relative velovity = ?)

%% User defined data - LEO Satellite parameters
% Satellite altitude 552 [km]
% Satellite (EIRP) = 34 [dBW/MHz]
% Satellite antenna gain 30 [dBi]
% Equivalent satellite antenna aperture 2 [m]
% Antenna gain-to-noise-temperature (G/T) 1.1 [dB/K]
SAT_EIRP = 34;  % [dBW/MHz]
% EIRP_density = calculate_EIRP_density(EIRP, bandwidth_kHz);
fprintf("EIRP Density = %f (dBW/MHz)\n", SAT_EIRP)
G_T = 1.1;

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
% satellite : omnidirectional
% since a parabolic antenna has a large volume, research on beamforming through a phased array antenna is being conducted
% ground stations : phased array antenna 相位陣列天線(iridium : custom 48-beam)

% **Regulatory Considerations** 
% S-band (2 to 4 GHz) 
% commonly used for satellite communication has allocated frequencies for such purposes by regulatory bodies like ITU
% Researchers may choose this band because of its availability and regulatory approval for satellite communication systems

% fq represents the frequency of the communication signal in Hertz (Hz)
fc = 2e9;       % center frequency in Hz(2GHz -> Hz)
fc_GHz = fc / 1e9;
fprintf("center frequency = %d GHz\n", fc_GHz)

% txpower represents the transmit power of the transmitter
% dBW (dBm=dBW+30)
% txpower = 200mW = dBm = dBW
txpower = 200; %[mW]
fprintf("transmit power = %d mW\n", txpower)

% beam_width represents the beam width of the phased array antenna in degrees
% This parameter determines the angular width of the main lobe of each beam
% 相控陣天線的波束寬度（degree）決定每個波束主瓣的角寬度
beam_width = 5; % degrees

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
fprintf("EIRP = %f (dBW)\n", EIRP) % -6.99

%% User defined data - Ground stations locaiton
Ground_Station();
% inclination (elevation angle)發射傾角 degrees (>10 degree)
alpha = inclination;       
% phi = calculate_phi(gs_lat, gs_lon); 
% central angle (地心_天北極 & 地心與衛星連線夾角)
% distance between sat & UE (meters)
d = calculate_d(Re, alpha, h0);
fprintf("communication path distance = %f (km)\n", d/1000)

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
% Boltzmann's constant [k] −228.6 dBW/K/Hz
% Bandwidth [B] 52.55 dBHz
bandwidth_dBHz = 52.55; %[dBHz]
bandwidth_kHz = 180;    %[kHz]

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
% Bandwidth = 750e3; % Bandwidth in Hz
BitRate = 512e3; % Bit Rate in Hz (bit/s)

%% Calculte free-space path loss
FSPL = FSPL_model(fc_GHz, d);
fprintf("FSPL = %f (dB)\n", FSPL)

%% Calculate path loss
% PL = path loss
% FSPL = free space path loss
% APL = atmospheric loss
% SMPL = shadowing margin loss
% SLPL = scintillation loss
% ADPL = additional loss
% PL = FSPL + APL + SMPL + SLPL + ADPL
PL = path_loss_model(FSPL);

%% Doppler shift model
fq_ref = fc;  % f reference
fq_edge = fc; % f edge
fq_res = fc;  % f residual
fd = Doppler_model(fc, v, alpha);
fprintf("Doppler model : %f (kHz)\n", fd)

%% transmitter & receiver antenna (SAT <-> GS)
% transmitter & receiver antenna(starlink) 
% GS : Parabolic antenna (dish antenna)
% SAT : phased array antenna

% antenna_transmission(sc, starlink, fc, txpower, PL)

% antennaType = "omnidirectional";
% antennaType = "Custom 48-Beam";

% if antennaType == "Custom 48-Beam" 
%   antenna = helperCustom48BeamAntenna(fq);
%   tx = transmitter(sats, ...
%        Frequency=fq, ...
%        Power=txpower, ...
%        Antenna=antenna); 
%        isotropic = arrayConfig(Size=[1 1]);
%        isotropicAntenna(tx,Array= isotropic);
% end
% isotropic = arrayConfig(Size=[1 1]);
% rx = receiver(gs,Antenna=isotropic);

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
SNR_spec = EbNo_To_SNR(EbNo_spec, BitRate, bandwidth_dBHz) ;

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
% SCS = 15;
% Multiple access = SC-FDMA
% Modulation QPSK
% Channel coding 1/3 Turbo code
TBS = 120;    %(bits)
% Fading channel model 3GPP TDL-D
% Residual Doppler shift 0, 950 [Hz]
NRep = 2;

%% 3.2 Link performance
% HARQ Configuration: 
% Configure parameters related to Hybrid Automatic Repeat reQuest (HARQ) 
% such as the number of retransmissions, coding scheme, redundancy versions 
% HARQ plays a crucial role in improving link reliability and performance

% Define Communication Scenarios: Define various communication scenarios considering different numbers of LEO satellites
% communication time variations, RTT (Round-Trip Time), packet sizes, and data transmission rates
% This will help in assessing the impact of satellite communication on link performance

%% Uplink & Downlink model
UL_model(SAT_EIRP, G_T, SAT_altitude)
DL_model(SAT_EIRP, G_T, SAT_altitude)

%% CNR calculate
% CNR = EIRP - FSPL + G_T - 228.6 - bandwidth
% EIRP : dBW, -7
% FSPL : dB, 153.31
% G_T : antenna-gain-to-noise-temperature [dB/k], 1.1
% Boltzman constant = 228.6 [dB]
% Bandwidth = 52.55 [dBHz], 180 [kHz]
data_rate = 0.01;        % 10 kbps = 0.01 Mbps
txpower_dBW = mW_to_dBW(txpower);
CNR = calculate_CNR(txpower_dBW ,fc_GHz, Transmit_antenna_gain, d, G_T, data_rate);

%% Ch4 Proposed solution
% 4.1 Reduction in  beam coverage
% 4.2 Addition of DMRS symbol

%% Ch5 Performance Analysis
% 5.1 Block Error Rate (BLER)
% 5.2 Throughput
% 5.3 Link Margin

%% Relay Algorithm
% Relay_algorithm()

%% 自定義衛星
% LEO_Sat_1 = satellite(sc, semi_major_axis, eccentricity, inclination, raan, arg_of_periapsis, trueAnomaly, Name = "LEO 1", OrbitPropagator = "two-body-keplerian");
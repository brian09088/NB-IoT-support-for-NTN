
- Channel Model 
    - fading, pathloss, modulation, uplink, downlink
1.  **Basic setup**
    - Satellite altitude = 650 km
    - communication path distance = 1597.492286 (km)
    - EIRP Density = 34.000000 (dBW/MHz)
    - center frequency = 2 GHz
    - transmit power = 200 mW
    - UE_antenna_gain = 0
    - Sat_antenna_gain = 30
    - Transmit_antenna_gain = 0
    - UE_antenna_power = 200.000000 mW
    - UE_antenna_power = 52.983174 dBm
2.  **others** 
    - EIRP = -6.989700 (dBW)
    - FSPL = 222.539375 (dB)
    - Doppler model : -0.000000 (kHz) 
    - Fading channel model : 3GPP TDL-D
    - maximum possible Doppler shift compensation (fd_max) = 1 (kHz)
------
- Reference CNR: 9.3198 dB 
- Uplink

 | ElevationAngle(degrees) | CNR(dB)|  FSPL(dB) |
 | ----------------------- | ------- | -------- |
 |           5             |   -4.4943   | 166.24 |
 |          15             |   -1.4878   | 163.23 |
 |          30             |    2.0158   | 159.73 |
 |          60             |    5.9015   | 155.84 |
 |          90             |    7.0198   | 154.73 |
  
- Downlink

 | ElevationAngle(degrees) | Link Margin(dB)| NRep_Add |
 | ----------------------- | ------- | -------- |
 |           5             |   -13.814   | 24 |
 |          15             |   -10.808   | 12 |
 |          30             |    -7.304   | 5 |
 |          60             |   -3.4183   | 2 |
 |          90             |   -2.3   | 1 |
 
------ 
- Reference CNR: 0.2889 dB
- Uplink 

 | ElevationAngle(degrees) | CNR(dB)|  FSPL(dB) |
 | ----------------------- | ------- | -------- |
 |           5             |   -3.6656   | 166.24 |
 |          15             |   -0.65901  | 163.23 |
 |          30             |    2.8445   | 159.73 |
 |          60             |    6.7303   | 155.84 |
 |          90             |    7.8485   | 154.73 |
  
- Downlink

 | ElevationAngle(degrees) | Link Margin(dB)| NRep_Add |
 | ----------------------- | ------- | -------- |
 |           5             |   -3.9545   | 2 |
 |          15             |   -0.9479   | 1 |
 |          30             |   2.5556   | 0 |
 |          60             |   6.4414   | 0 |
 |          90             |   7.5596   | 0 |
 
------ 
 - 因為衛星天線屬於射頻晶片（RF），所以採用SNR而非SNR
 - SNR是總訊號功率比雜訊。 CNR有用訊號功率比噪音
 - **Satellite CNR Config**:
     - TransmitterPower: 22.9832
     - TransmitterSystemLoss: 0
     - TransmitterAntennaGain: 0
     - Distance: 1.5975e+03
     - Frequency: 2
     - MiscellaneousLoss: 0
     - GainToNoiseTemperatureRatio: 1.1000
     - ReceiverSystemLoss: 0
     - BitRate: 0.0100
     - SymbolRate: 10
     - Bandwidth: 6
- **CN = 22.3637 dB**
    -  TransmitterEIRP: 22.9832
    - FSPL: 162.5372
    - ReceivedIsotropicPower: -139.5540
    - CarrierToNoiseDensityRatio: 90.1452
    - ReceivedEbNo: 50.1452
    - ReceivedEsNo: 20.1452 

------

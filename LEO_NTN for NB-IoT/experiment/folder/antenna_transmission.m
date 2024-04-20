% function for showing map coverage and ground stations locaiton
% also contain antenna transmit and receive (satllite and ground stations)

function [] = antenna_transmission(sc, starlink, fc, txpower, PL)
    
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

    antennaType = "Uniform Rectangular Array";
    
    % Satellite transmitter :
    % starlink : phased array antenna
    if antennaType == "Uniform Rectangular Array"
        
        % Define parameters for the phased array antenna
        numElements = 64;  % Number of elements in the array
        elementSpacing = 0.5;  % Spacing between elements in wavelengths
        phasedArray = phased.URA('Size',[8,8],'ElementSpacing',[elementSpacing,elementSpacing]);
        
        % Define transmitter with the phased array antenna
        tx = transmitter(starlink, ...
            'Frequency', fc, ...
            'MountingAngles', [0,-90,0], ... % [yaw, pitch, roll] with -90 using Phased Array System Toolbox convention
            'Power', txpower, ...
            'Antenna', phasedArray);

    end
    
    % Add Ground Station Receivers (地面接收站獲取訊號後透過基地台轉傳)
    % GS : Parabolic antenna (dish antenna)
    % antennaType = "Gaussian Antenna";
    
    gimbalrxSat = gimbal(starlink);
    gimbaltxSat = gimbal(starlink);

    gainToNoiseTemperatureRatio = 5;                                                        % dB/K
    systemLoss = 3;                                                                         % dB
    rxSat = receiver(gimbalrxSat,Name="Satellite Receiver",GainToNoiseTemperatureRatio= ...
        gainToNoiseTemperatureRatio,SystemLoss=systemLoss);

    % s-band 2GHz
    % txpower = 200mW = 23 dBm = -7 dBW
    % bitRate ranging from 50 Mbps to 150 Mbps
    frequency = 2e9;                                                                     % Hz
    power = -7;                                                                            % dBW
    bitRate = 100;                                                                         % Mbps
    systemLoss = PL;                                                                       % dB                                                                     % dB
    txSat = transmitter(gimbaltxSat,Name="Satellite Transmitter",Frequency=frequency, ...
        power=power,BitRate=bitRate,SystemLoss=systemLoss);
    
    % 594mm(Starlink Dish diameter)
    dishDiameter = 0.5;                                                                    % meters
    apertureEfficiency = 0.5;
    gaussianAntenna(txSat,DishDiameter=dishDiameter,ApertureEfficiency=apertureEfficiency);
    gaussianAntenna(rxSat,DishDiameter=dishDiameter,ApertureEfficiency=apertureEfficiency);

    gs1 = groundStation(sc, 25.02, 121.31, Name="GS 1");
    gs2 = groundStation(sc, 22.36, 120.18, Name="GS 2");

    pointAt(gimbaltxSat,gs2);
    pointAt(gimbalrxSat,gs1);

    gimbalgs1 = gimbal(gs1);
    gimbalgs2 = gimbal(gs2);
    
    frequency = 2e9;                                                                          % Hz
    power = -7;                                                                                % dBW
    bitRate = 100;                                                                              % Mbps
    txGs1 = transmitter(gimbalgs1,Name="GS 1 Transmitter",Frequency=frequency, ...
            Power=power,BitRate=bitRate);

    requiredEbNo = 14;                                                                     % dB
    rxGs2 = receiver(gimbalgs2,Name="Gs 2 Receiver",RequiredEbNo=requiredEbNo);
    
    dishDiameter = 0.5;                                % meters
    gaussianAntenna(txGs1,DishDiameter=dishDiameter);
    gaussianAntenna(rxGs2,DishDiameter=dishDiameter);
    
    pointAt(gimbalgs1,starlink);
    pointAt(gimbalgs2,starlink);

    lnk = link(txGs1,rxSat,txSat,rxGs2);
    % antenna_interference();
    
    % ground station coverage/distance/numbers...
    % have some problem to be fix
    % pattern(tx,Size=500000);

end
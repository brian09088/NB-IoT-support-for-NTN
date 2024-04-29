clc,clear,close all
% Create and Visualize Scenario
startTime = datetime(2024,2,26,12,0,0);
stopTime = startTime + minutes(30);
sampleTime = 5; % seconds
sc = satelliteScenario(startTime,stopTime,sampleTime);
viewer = satelliteScenarioViewer(sc,ShowDetails=false);

% Add Satellite Constellation
% As of February 2024, there are 5,438 Starlink satellites in orbit
numSatellitesPerOrbitalPlane = 11;
numOrbits = 6;

RAAN = [];
trueanomaly = [];

for i = 1:numOrbits
    for j = 1:numSatellitesPerOrbitalPlane 
        
        RAAN(end+1) = 180*(i-1)/numOrbits; %#ok<SAGROW>
        if mod(i,2)
            trueanomaly(end+1) = 360*(j-1)/numSatellitesPerOrbitalPlane; %#ok<SAGROW>
        else
            % Satellites offset in alternating orbits
            trueanomaly(end+1) = 360*(j-1 + 0.5)/numSatellitesPerOrbitalPlane; %#ok<SAGROW> 
        end
    end
end

% 衛星基本參數配置
semimajoraxis = repmat((6378 + 780)*1e3,size(RAAN)); % 低軌衛星高度780, kms
inclination = repmat(86.4,size(RAAN)); % degrees
eccentricity = zeros(size(RAAN)); % degrees
argofperiapsis = zeros(size(RAAN)); % degrees

sats = satellite(sc,semimajoraxis,eccentricity,inclination,RAAN,argofperiapsis,trueanomaly,Name="LEO_NTN_" + string(1:66)');

% Add Grid of Ground Stations Covering Taiwan-ROC(地面基地台)
latlim = [22.00417 25.12825];
lonlim = [118.31833 121.753];
spacingInLatLon = 1; % 基站間距離，越小覆蓋率越高但計算時間較長(單位degrees)
EPSG_code = 3826 % (TWD97-EPSG:3826/3830)
tw_lat = [20 27.5];
tw_lon = [115 125];
proj = projcrs(EPSG_code); % projected coordinate reference system
spacingInXY = deg2km(spacingInLatLon)*1000; % 相對於spacingInLatLon轉換為xy座標系統(單位:meters)
[xlim,ylim] = projfwd(proj,latlim,lonlim);
R = maprefpostings(xlim,ylim,spacingInXY,spacingInXY);
[X,Y] = worldGrid(R);
[gridlat,gridlon] = projinv(proj,X,Y);

% 以下為copilot提供 GCS座標系統->PCS座標系統 轉換方法:

% 您看到的錯誤訊息是因為 MATLAB 中的 projcrs 函數需要投影座標參考系統 (CRS) 程式碼作為參數1。 您使用的 EPSG 代碼 (3826) 適用於地理座標系，而非投影座標系。
% 投影 CRS 提供將笛卡爾 x 和 y 地圖座標分配給物理位置的資訊。 投影 CRS 由地理 CRS 和幾個用於在地理 CRS1 與座標之間進行座標轉換的參數組成。
% 您可以透過多種方式建立投影 CRS 物件1：
% 使用readgeoraster匯入柵格數據，然後查詢傳回的柵格參考對象的ProjectedCRS屬性。
% 使用readgeotable函數匯入向量數據，然後在傳回的地理空間表中查詢形狀物件的ProjectedCRS屬性。
% 使用 shapeinfo 函數取得有關 shapefile 的信息，然後查詢傳回的結構體的 CooperativeReferenceSystem 欄位。
% 使用 projcrs 函數（此處描述）。
% 您需要為您的特定用例找到正確的預計 CRS 代碼。 如果您不確定代碼是否正確，您可能需要諮詢 EPSG 註冊表或其他相關資源。
% 如果您仍然遇到問題，您可能需要考慮使用 try/catch 語句來處理錯誤2。 這可以幫助您識別和管理程式碼中的潛在問題。
% 請注意，此建議是基於 MATLAB 文件中提供的資訊以及您提供的程式碼。 對於更具體或更複雜的問題，您可能需要諮詢 MathWorks 支援或 MATLAB 使用者社群。

% Define the EPSG code for TWD97
% EPSG_code = 3826;
% Create a projected CRS object using the specified EPSG code
% proj = projcrs(EPSG_code);
% Now you can use the 'proj' object for further transformations or calculations
% related to the TWD97 coordinate system.
% Example: Unproject x-y map coordinates to latitude-longitude coordinates
% (assuming you have x and y coordinates in TWD97)
% Transform the map coordinates to latitude-longitude coordinates
% [lat, lon] = projinv(proj, x, y);
% Proceed with the rest of your simulation setup using the 'proj' object and
% the transformed latitude-longitude coordinates.
% Note: Replace '...' with your actual x and y coordinates.
% Rest of your simulation code...
% (Add the remaining part of your satellite constellation setup and other components)

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

%% 設計所需的天線陣列與波束
% satellite : phased array antenna (SPAA)
% gs : phased array antenna (see GS_antenna.m)
fq = 1e9; % Hz
txpower = 30; % dBW
antennaType = "Custom 48-Beam";
beamWidth = 5; % degrees

% Custom-48 Beam 天線：
% Custom-48 Beam 天線是一種具有 48 個波束的天線。
% 這種天線通常用於 LEO 衛星的通信，例如像 Iridium NEXT 衛星這樣的衛星。
% 每個波束可以獨立調整，以實現不同方向的波束形成。
% 在 MATLAB 中，您可以使用 Phased Array System Toolbox 來模擬這種天線。
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

% Visualize Coverage on an axesm-Based Map
minpowerlevel = -120; % dBm
maxpowerlevel = -100; % dBm

figure
worldmap(latlim,lonlim)
mlabel south

colormap turbo
clim([minpowerlevel maxpowerlevel])
geoshow(gridlat,gridlon,maxsigstrength,DisplayType="contour",Fill="on")
geoshow(tw_lat,tw_lon,Facecolor="none")
% taiwanBoundaries = geoshape([...
%    20.5170, 118.2008; % Southwest corner
%    25.3000, 122.0000; % Northwest corner
%    25.3000, 121.3000; % Northeast corner
%    20.6834, 119.6000; % Southeast corner
%    20.5170, 118.2008]);  % Repeat the first point to close the polygon
% geoshow(taiwanBoundaries, FaceColor = "none");

cBar = contourcbar;
title(cBar,"dBm");
title("Signal Strength at " + string(startTime) + " UTC")

% 設定當地重要城市經緯度，進行定位與標註
TP_location = [25.02 121.31]; % Taipei
TC_location = [24.08 120.41]; % Taichung
Kao_location = [22.36 120.18]; % Kaosiung

textm(TP_location(1),TP_location(2),"Taipei")
textm(TC_location(1),TC_location(2),"Taichung")
textm(Kao_location(1),Kao_location(2),"Kaohsiung")

% Compute Coverage Map Contours (NB-IoT Coverage Enhancement levels 部分待調整...)
levels = linspace(minpowerlevel,maxpowerlevel,11);
GT = contourDataGrid(gridlat,gridlon,maxsigstrength,levels,proj);
GT = sortrows(GT,"Power (dBm)");
disp(GT)

% Visualize Coverage on a Map Axes
figure
newmap(proj)
hold on

colormap turbo
clim([minpowerlevel maxpowerlevel])
geoplot(GT,ColorVariable="Power (dBm)",EdgeColor="none")
% taiwanBoundaries = geoshape([...
%     20.5170, 118.2008; % Southwest corner
%     25.3000, 122.0000; % Northwest corner
%     25.3000, 121.3000; % Northeast corner
%     20.6834, 119.6000; % Southeast corner
%     20.5170, 118.2008  % Repeat the first point to close the polygon]);
geoplot(tw_lat,tw_lon, 'FaceColor', 'none', 'EdgeColor', 'b');

cBar = colorbar;
title(cBar,"dBm");
title("Signal Strength at " + string(startTime) + " UTC")

text(TP_location(1),TP_location(2),"Taipei",HorizontalAlignment="center")
text(TC_location(1),TC_location(2),"Taichung",HorizontalAlignment="center")
text(Kao_location(1),Kao_location(2),"Kaohsiung",HorizontalAlignment="center")

% Compute and Visualize Coverage for a Different Time
secondTOI = startTime + minutes(2); % 2 minutes after the start of the scenario
maxsigstrength = satcoverage(gridpts,sc,secondTOI,inregion,beamWidth);

GT2 = contourDataGrid(gridlat,gridlon,maxsigstrength,levels,proj);
GT2 = sortrows(GT2,"Power (dBm)");

figure
newmap(proj)
hold on

colormap turbo
clim([minpowerlevel maxpowerlevel])
geoplot(GT2,ColorVariable="Power (dBm)",EdgeColor="none")
% taiwanBoundaries = geoshape([...
%     20.5170, 118.2008; % Southwest corner
%     25.3000, 122.0000; % Northwest corner
%     25.3000, 121.3000; % Northeast corner
%     20.6834, 119.6000; % Southeast corner
%     20.5170, 118.2008  % Repeat the first point to close the polygon]);
geoplot(tw_lat,tw_lon, 'FaceColor', 'none', 'EdgeColor', 'b');

cBar = colorbar;
title(cBar,"dBm");
title("Signal Strength at " + string(secondTOI) + " UTC")

text(TP_location(1),TP_location(2),"Taipei",HorizontalAlignment="center")
text(TC_location(1),TC_location(2),"Taichung",HorizontalAlignment="center")
text(Kao_location(1),Kao_location(2),"Kaohsiung",HorizontalAlignment="center")

% Compute Coverage Levels for Cities
covlevels1 = [coveragelevel(sydlatlon(1),sydlatlon(2),GT); ...
    coveragelevel(mellatlot(1),mellatlot(2),GT); ...
    coveragelevel(brislatlon(1),brislatlon(2),GT)];
covlevels2 = [coveragelevel(sydlatlon(1),sydlatlon(2),GT2); ...
    coveragelevel(mellatlot(1),mellatlot(2),GT2); ...
    coveragelevel(brislatlon(1),brislatlon(2),GT2)];

covlevels = table(covlevels1,covlevels2, ...
    RowNames=["Taipei" "Taichung" "Kaohsiung"], ...
    VariableNames=["Signal Strength at T1 (dBm)" "Signal Strength T2 (dBm)"]);

% Helper Functions
function coveragedata = satcoverage(gridpts,sc,timeIn,inregion,beamWidth)

    % Get satellites and ground station receivers
    sats = sc.Satellites;
    rxs = [sc.GroundStations.Receivers];

    % Compute the latitude, longitude, and altitude of all satellites at the input time
    lla = states(sats,timeIn,"CoordinateFrame","geographic");

    % Initialize coverage data
    coveragedata = NaN(size(gridpts));

    for satind = 1:numel(sats)
        % Create a geopolyshape for the satellite field-of-view
        fov = fieldOfViewShape(lla(:,1,satind),beamWidth);

        % Find grid and rx locations which are within the field-of-view
        gridInFOV = isinterior(fov,gridpts);
        rxInFOV = gridInFOV(inregion);

        % Compute sigstrength at grid locations using temporary link objects
        gridsigstrength = NaN(size(gridpts));
        if any(rxInFOV)
            tx = sats(satind).Transmitters;
            lnks = link(tx,rxs(rxInFOV));
            rxsigstrength = sigstrength(lnks,timeIn)+30; % Convert from dBW to dBm
            gridsigstrength(inregion & gridInFOV) = rxsigstrength;
            delete(lnks)
        end

        % Update coverage data with maximum signal strength found
        coveragedata = max(gridsigstrength,coveragedata);
    end
end

function satFOV = fieldOfViewShape(lla,beamViewAngle)

    % Find the Earth central angle using the beam view angle
    if isreal(acosd(sind(beamViewAngle)*(lla(3)+earthRadius)/earthRadius))
        % Case when Earth FOV is bigger than antenna FOV 
        earthCentralAngle = 90-acosd(sind(beamViewAngle)*(lla(3)+earthRadius)/earthRadius)-beamViewAngle;
    else
        % Case when antenna FOV is bigger than Earth FOV 
        earthCentralAngle = 90-beamViewAngle;
    end

    % Create a buffer zone centered at the position of the satellite with a buffer of width equaling the Earth central angle
    [latBuff,lonBuff] = bufferm(lla(1),lla(2),earthCentralAngle,"outPlusInterior",100);

    % Handle the buffer zone crossing over -180/180 degrees
    if sum(abs(lonBuff) == 180) > 2
        crossVal = find(abs(lonBuff)==180) + 1;
        lonBuff(crossVal(2):end) = lonBuff(crossVal(2):end) - 360 *sign(lonBuff(crossVal(2)));
    elseif sum(abs(lonBuff) == 180) == 2
        lonBuff = [lonBuff; lonBuff(end); lonBuff(1); lonBuff(1)];
        if latBuff(1) > 0
            latBuff = [latBuff; 90; 90; latBuff(1)];
        else
            latBuff = [latBuff; -90; -90; latBuff(1)];
        end
    end

    % Create geopolyshape from the resulting latitude and longitude buffer zone values
    satFOV = geopolyshape(latBuff,lonBuff);
end

function GT = contourDataGrid(latd,lond,data,levels,proj)

    % Pad each side of the grid to ensure no contours extend past the grid bounds
    lond = [2*lond(1,:)-lond(2,:); lond; 2*lond(end,:)-lond(end-1,:)];
    lond = [2*lond(:,1)-lond(:,2) lond 2*lond(:,end)-lond(:,end-1)];
    latd = [2*latd(1,:)-latd(2,:); latd; 2*latd(end,:)-latd(end-1,:)];
    latd = [2*latd(:,1)-latd(:,2) latd 2*latd(:,end)-latd(:,end-1)];
    data = [ones(size(data,1)+2,1)*NaN [ones(1,size(data,2))*NaN; data; ones(1,size(data,2))*NaN] ones(size(data,1)+2,1)*NaN];

    % Replace NaN values in power grid with a large negative number
    data(isnan(data)) = min(levels) - 1000;
    
    % Project the coordinates using an equal-area projection
    [xd,yd] = projfwd(proj,latd,lond);
    
    % Contour the projected data
    fig = figure("Visible","off");
    obj = onCleanup(@()close(fig));
    c = contourf(xd,yd,data,LevelList=levels);
    
    % Initialize variables
    x = c(1,:);
    y = c(2,:);
    n = length(y);
    k = 1;
    index = 1;
    levels = zeros(n,1);
    latc = cell(n,1);
    lonc = cell(n,1);
    
    % Calculate the area within each contour line. Remove areas that
    % correspond to holes and ignore negative areas.
    while k < n
        level = x(k);
        numVertices = y(k);
        yk = y(k+1:k+numVertices);
        xk = x(k+1:k+numVertices);
        k = k + numVertices + 1;
    
        [first,last] = findNanDelimitedParts(xk);
        jindex = 0;
        jy = {};
        jx = {};
        sumpart = 0;
    
        for j = 1:numel(first)
            % Process the j-th part of the k-th level
            s = first(j);
            e = last(j);
            cx = xk(s:e);
            cy = yk(s:e);
            if cx(end) ~= cx(1) && cy(end) ~= cy(1)
                cy(end+1) = cy(1); 
                cx(end+1) = cx(1);
            end

            if j == 1 && ~ispolycw(cx,cy)
                % First region must always be clockwise
                [cx,cy] = poly2cw(cx,cy);
            end
    
            jindex = jindex + 1;
            jy{jindex,1} = cy(:)';
            jx{jindex,1} = cx(:)';
    
            a = polyarea(cx,cy);
            if ~ispolycw(cx,cy)
                % Remove areas that correspond to holes
                a = -a;
            end

            sumpart = sumpart + a;
        end
    
        % Add a part when its area is greater than 0. Unproject the
        % coordinates.
        [jx,jy] = polyjoin(jx,jy);
        if length(jy) > 2 && length(jx) > 2 && sumpart > 0
            [jlat,jlon] = projinv(proj,jx,jy);
            latc{index,1} = jlat(:)';
            lonc{index,1} = jlon(:)';
            levels(index,1) = level;
            index = index + 1;
        end
    end
    
    % Create contour shapes from the unprojected coordinates. Include the
    % contour shapes, the areas, and the power levels in a geospatial
    % table.
    latc = latc(1:index-1);
    lonc = lonc(1:index-1);
    Shape = geopolyshape(latc,lonc);
    Levels = levels(1:index-1);
    allPartsGT = table(Shape,Levels);  

    % Condense the geospatial table into a new geospatial table with one
    % row per power level.
    GT = table.empty;
    levels = unique(allPartsGT.Levels);
    for k = 1:length(levels)
        gtLevel = allPartsGT(allPartsGT.Levels == levels(k),:);
        tLevel = geotable2table(gtLevel,["Latitude","Longitude"]);
        [lon,lat] = polyjoin(tLevel.Longitude,tLevel.Latitude);
        Shape = geopolyshape(lat,lon);
        Levels = levels(k);
        GT(end+1,:) = table(Shape,Levels);
    end

    maxLevelDiff = max(abs(diff(GT.Levels)));
    LevelRange = [GT.Levels GT.Levels + maxLevelDiff];
    GT.LevelRange = LevelRange;

    % Clean up the geospatial table
    GT.Properties.VariableNames = ...
        ["Shape","Power (dBm)","Power Range (dBm)"]; 
end

function powerLevels = coveragelevel(lat,lon,GT)

    % Determine whether points are within coverage
    inContour = false(length(GT.Shape),1);
    for k = 1:length(GT.Shape)
        inContour(k) = isinterior(GT.Shape(k),geopointshape(lat,lon));
    end

    % Find the highest coverage level containing the point
    powerLevels = max(GT.("Power (dBm)")(inContour));

    % Return -inf if the point is not contained within any coverage level
    if isempty(powerLevels)
        powerLevels = -inf;
    end
end

function [first,last] = findNanDelimitedParts(x)

    % Find indices of the first and last elements of each part in x. 
    % x can contain runs of multiple NaNs, and x can start or end with 
    % one or more NaNs.

    n = isnan(x(:));
    
    firstOrPrecededByNaN = [true; n(1:end-1)];
    first = find(~n & firstOrPrecededByNaN);
    
    lastOrFollowedByNaN = [n(2:end); true];
    last = find(~n & lastOrFollowedByNaN);
end
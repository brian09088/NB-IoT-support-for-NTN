%% function for convert satellite TLE to coordinate

function [location] = satelliteTLE(sat_csv)
    
    startTime = datetime(2021,11,11,8,30,0);
    sampleTime = 30;              % seconds
    stopTime = startTime + sampleTime;      
    sc = satelliteScenario(startTime,stopTime,sampleTime);

    sat = satellite(sc,sat_csv);
    time = datetime(2021,11,11,9,30,0);
    location = states(sat,time,"CoordinateFrame","geographic");
    
end
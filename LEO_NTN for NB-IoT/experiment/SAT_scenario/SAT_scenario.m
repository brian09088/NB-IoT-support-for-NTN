classdef SAT_scenario
    
    properties
        % get from main.m 
        sc
        starlink
        sat
        fc
        txpower
        PL
        % antenna
        antennaType
        numElements
        elementSpacing
        % for scenario observation
        gainToNoiseTemperatureRatio
        systemLoss
        frequency  
        power               
        bitRate                    
        dishDiameter
        apertureEfficiency
        requiredEbNo
    end
    
    methods
        % Construct
        function obj = UE(id)
            % Construct an instance of this class
            obj.gainToNoiseTemperatureRatio = 5;
            obj.systemLoss = 3;
            obj.frequency = 2e9; 
            obj.power = -7;              
            obj.bitRate = 100;                   
            obj.dishDiameter = 5;     % meters
            obj.apertureEfficiency =0.5;
            obj.requiredEbNo = 14;
        end
    end
end


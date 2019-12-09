classdef cgmObj
    % This object represents the continuous glucose monitor. The aim is to
    % take real-time values of the plasma glucose conentration and the
    % plasma insulin cocncentration. This accounts for the parametric
    % delays in reporting by adding a noise signal to the model.
    
    properties
        sensor_MeasurementDelay {mustBeNumeric}
    end
    
    methods
        function obj = untitled2(inputArg1,inputArg2)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end


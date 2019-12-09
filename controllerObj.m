classdef controllerObj
    % This is the basic template for creating controllers for handling the
    % control problem at hand, "T1DM".
    
    properties
        controller_Status
        optimal_dose
    end
    
    methods
        function obj = untitled(inputArg1,inputArg2)
            % This displays a message about the type of controller created
            % and its status
            
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


classdef calibrateDLRs
    
    properties
        Property1
    end
    
    methods
        function self = calibrateDLRs(inputArg1)
            self.Property1 = inputArg1;
        end
        
        function outputArg = method1(obj,inputArg)
            outputArg = obj.Property1 + inputArg;
        end
    end
end


classdef Config
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function res = defaultConfig()
            %  templateFile = name of template file, in dataPathIn
            %               [defaults to 'template.txt']
            %  basisSet = basis set keyword (Gaussian format)
            %               [defaults to 'STO-3G']
            %  method = Method that you want to use
            %              [defults to hf]
            %  charge = charge on the fragment
            %             [defaults to 0]
            %  spin   = spin (multiplicity) of the fragment,
            %             using Gaussian convention
            %             [defaults to 1]
            res.template = 'template';
            res.basisSet = 'STO-3G';
            res.method   = 'hf';
            res.charge   = 0;
            res.spin     = 1;
            res.calcEn   = 1;
            res.timeOut  = -1;
            res.zmat = ZMatrix();
        end
    end
    
end


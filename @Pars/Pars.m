classdef Pars < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %Where the par values are stored in the form of a list
        bond_pars = {};
        ang_pars = {};
        di_pars = {};
        
        bond_vals = {
            6 6 1.54
            6 1 1.14
            6 8 1.43
            6 7 1.47
            6 9 1.35
            6 17 1.77
            7 1 1.001
            7 8 1.25
            8 1 0.958
            };
        ang_vals = {
            120.00001   %Bent
            120.00001   %Tri
            109.5 %Tetrahedral
            };
        di_vals = {0.0001, 120.00001, -120.00001, 180.000001};
        
    end
    
    methods
    end
    
end


classdef Fitme2
    %FITME2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        parFitme
    end
    
    methods
        function res = err(obj,par)
            errs = cell(1, length(obj.parFitme));
            array = obj.parFitme;
            parfor i = 1:length(obj.parFitme)
                errs{i} = array(i).err(par);
            end
            res = [];
            for i = 1:length(obj.parFitme)
                res = [res errs{i}];
            end
        end
        function res = getPars(obj)
            res = zeros(1,obj.parFitme(1).npar);
            ic = 1;
            for i = 1:size(obj.parFitme(1).mixers,2)
                mtemp = obj.parFitme(1).mixers{1,i};
                np = mtemp.npar;
                if (np > 0)
                    res(ic:(ic+np-1)) = mtemp.getPars;
                end
                ic = ic + np;
            end
        end
        function setPars(obj,par)
            % sets parameters, and updates densities
            if (size(par,2) ~= obj.npar)
                error(['Fitme.SetPars called with ',num2str(size(par,2)), ...
                    ' parameters when ',num2str(obj.npar),' are needed ']);
            end
            for j = 1:length(obj.parFitme)
                ic = 1;
                mix = obj.parFitme(j).mixers;
                for i = 1:size(mix,2)
                    mtemp = mix{1,i};
                    np = mtemp.npar;
                    if (np > 0)
                        mtemp.setPars( par(ic:(ic+np-1)));
                    end
                    ic = ic + np;
                end
            end
        end
    end
    
end


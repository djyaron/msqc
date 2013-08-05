function runFragment( fragment )

            nparIn = size(fragment.config.par,1) * size(fragment.config.par,2);
            if (nparIn ~= fragment.npar)
                error(['template has ',num2str(fragment.npar),' parameters',...
                    ' while config contains ',num2str(nparIn),' pars']);
            end

            %%
            [found,fragment.fileprefix] = Fragment.findCalc(fragment.dataPath,fragment.config);
            if (~found)
                createTempFolder( fragment );
                save([fragment.fileprefix,'_cfg.mat'], 'fragment.config' );
                fragment.initializeZipData( [fragment.fileprefix,'.zip'] );
            end
            fragment.loadZipData( [fragment.fileprefix,'.zip'] );
            
            fragment.nenv = 0;
            % Set the environment array to have the correct class type
            if exist( 'Environment', 'class' ) == 8
                fragment.env = Environment.empty(0,0);
            end

end

function createTempFolder( fragment )
    temp1 = tempname('a'); % makes "a\uniquestring"
    uniqueStr = temp1(3:end);
    fragment.fileprefix = [fragment.dataPath,filesep,fragment.config.template, ...
        '_',uniqueStr];
end
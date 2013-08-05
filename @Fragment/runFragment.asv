function runFragment( fragment, dataPathIn, configIn )

    if (nargin < 2)
        fragment.dataPath = 'data';
        dataPath = 'data';
    else
        fragment.dataPath = dataPathIn;
    end
    if (nargin < 3)
        config = Config();
        fragment.config = config;
    else
        fragment.config = configIn;
    end
    
    nparIn = size(fragment.config.par,1) * size(fragment.config.par,2);
    if (nparIn ~= fragment.npar)
        error(['template has ',num2str(fragment.npar),' parameters',...
            ' while config contains ',num2str(nparIn),' pars']);
    end

    %%
    [found,fragment.fileprefix] = Fragment.findCalc(dataPath,config);
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
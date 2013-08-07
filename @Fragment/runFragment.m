function bool = runFragment( fragment, dataPathIn, configIn )

    if (nargin < 2)
        fragment.dataPath = 'data';
        dataPath = 'data';
    else
        fragment.dataPath = dataPathIn;
        dataPath = dataPathIn;
    end
    if (nargin < 3)
        config = Config();
        fragment.config = configIn;
    else
        fragment.config = configIn;
        config = configIn;
    end

    %%
    [found,fragment.fileprefix] = Fragment.findCalc(dataPath,config);
    if (~found)
        createTempFolder( fragment );
        save([fragment.fileprefix,'_cfg.mat'], 'config' );
        bool = fragment.initializeZipData( [fragment.fileprefix,'.zip'] );
        if bool == 0
            return
        end
    end
    fragment.loadZipData( [fragment.fileprefix,'.zip'] );

    fragment.nenv = 0;
    % Set the environment array to have the correct class type
    if exist( 'Environment', 'class' ) == 8
        fragment.env = Environment.empty(0,0);
    end
    bool = 1;
end

function createTempFolder( fragment )
    temp1 = tempname('a'); % makes "a\uniquestring"
    uniqueStr = temp1(3:end);
    fragment.fileprefix = [fragment.dataPath,filesep,fragment.config.title, ...
        '_',uniqueStr];
end
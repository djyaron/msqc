function runFragment( fragment )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
            %  function comments
            %  dataPath = directory (including c:\ etc) for data storage
            %               do not include a \ at end of paths
            %               [defaults to 'data']
            %  configIn = configuration structure
            %               [defaults to 'Fragment.defaultConfig();
                       
%%          Num of arguments handling

            %Needs to be fixed cause of arg changes!!
            if (nargin < 1)
                fragment.dataPath = 'data';
            else
                fragment.dataPath = dataPathIn;
            end
            if (nargin < 2)
                fragment.config = Fragment.defaultConfig();
            else
                fragment.config = configIn;
            end
%%
            [found,fragment.fileprefix] = Fragment.findCalc(fragment.dataPath,fragment.config);
            % Backword compatibility (load mat file)
            
            %%              Pull config vars
            
            
            charge   = fragment.config.charge;
            spin     = fragment.config.spin;
            %%              Error Checks
            nparIn = size(fragment.config.par,1) * size(fragment.config.par,2);
            if (nparIn ~= fragment.npar)
                error(['template has ',num2str(fragment.npar),' parameters',...
                    ' while config contains ',num2str(nparIn),' pars']);
            end
            %%              Build GJF
            
            
            %%
            
            if (~found)
                temp1 = tempname('a'); % makes "a\uniquestring"
                uniqueStr = temp1(3:end);
                fragment.fileprefix = [fragment.dataPath,filesep,fragment.config.template, ...
                    '_',uniqueStr];
                % save config file, in *.mat format
                Cfile = fragment.config;
                save([fragment.fileprefix,'_cfg.mat'],  'Cfile' );
                % create and save zip file
                zipFile = [fragment.fileprefix,'.zip'];
                fragment.initializeZipData(zipFile);
            end
            zipFile = [fragment.fileprefix,'.zip'];
            fragment.loadZipData(zipFile);
            fragment.nenv = 0;
            % Set the environment array to have the correct class type
            if exist( 'Environment', 'class' ) == 8
                fragment.env = Environment.empty(0,0);
            end

end

function buildGjf( fragment )

    basisSet = fragment.config.basisSet;
    method   = fragment.config.method;
    
    headerObj = Header( basisSet, method );
    headerObj.link0 = {'rwf=temp.rwf' 'nosave' 'chk=temp.chk'}';
    if obj.config.opt == 1
        headerObj.route = {'opt'};
    end
    headerObj.output = {'nosymm int=noraff iop(99/6=1)' ...
        'scf=conventional' 'symm=noint'};
    fragment.config.header = headerObj;

    headerText = headerObj.makeHeader();
    title = [fragment.config.template, '\n\n'];
    charge_mult = [num2str(charge), ' ', num2str(spin), '\n'];
    zmat_body = zmat.build_gjf();

    fragment.gaussianFile = [headerText, title, charge_mult, zmat_body];
end
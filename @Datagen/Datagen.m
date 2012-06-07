classdef Datagen
    %DATAGEN Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % These vars must always be created first.
        jobname;        % Name of the data set. Ex. 'methane'.
        pars;           % Full cell array of parameters.
        
        % The remaining vars have defaults, but may be changed by taking
        % advantage of varargin in the apropriate function call.
        % 'root' implies a full path. 'path' may be full or relative to
        % msqcRoot.
        
        % Template Vars
        tplName;        % Filename of the main template file.
        tplNameGen;     % Filename of the gen template file.
        tplPath;        % Path to the main template file.
        tplPathGen;     % Path to the gen template file.
        
        % Path Vars
        msqcRoot;       % Full path to msqc code.
        dataPath;       % Path to data location.
        dataFolder;     % Name of folder inside dataPath to use.
        dataRoot;       % Full path to data.
        templatePath;   % Path to templates. It would be nice if this didn't 
                        % need to be in the data folder.
        
        % Environment Vars
        envType;        % 'charge' or 'field'.
        envFilename;    % Name of .mat file of envs.
        envFile;        % Path to .mat file.
        nenv;           % Number of environments.
        
        % Calculation Vars
        method;         % 'hf' or 'MP2'.
        HLbasis;        % Cell array of HL basis sets to use.
        datFile         % .mat file for HL, LL.
    end
    
    methods        
        function res = Datagen(jobname, pars, varargin)
            % Initialize a Datagen object with default or user specified
            % parameters.
            % jobname = Name of the data set. Ex. 'methane'.
            % pars = Full cell array of parameters.
            % varargin = specify nondefault values for any other class
            %            properties.
            
            % Set mandatory parameters.
            res.jobname = jobname;
            res.pars = pars;
            
            % Set options.
            res.msqcRoot = checkForInput(varargin, 'msqcRoot', ...
                'D:\Users\Alex\Programming\msqc\');
            res.dataPath = checkForInput(varargin, 'dataPath', ...
                ['data', filesep]);
            res.templatePath = checkForInput(varargin, 'templatePath', ...
                ['templates', filesep]);
            res.envFilename = checkForInput(varargin, 'envFilename', ...
                'env1.mat');
            res.tplName = checkForInput(varargin, 'tplName', ...
                jobname);
            res.tplNameGen = checkForInput(varargin, 'tplNameGen', ...
                [jobname, '-gen']);
            res.envType = checkForInput(varargin, 'envType', 'charge');
            res.nenv = checkForInput(varargin, 'nenv', 100);
            res.method = checkForInput(varargin, 'method', 'hf');
            res.HLbasis = checkForInput(varargin, 'HLbasis', {'6-31G'});
            
            % Set derivative vars.
            res.dataFolder = checkForInput(varargin, 'dataFolder', ...
                [res.dataPath, jobname, '_', res.envType, '_', res.method]);
            res.dataRoot = checkForInput(varargin, 'dataRoot', ...
                [res.msqcRoot, res.dataFolder]);
            res.datFile = checkForInput(varargin, 'datFile', ...
                [res.dataFolder, filesep, jobname, '_', res.envType, '_', ...
                res.method, '.mat']);
            res.tplPath = checkForInput(varargin, 'tplPath', ...
                [res.templatePath, res.tplName, '.tpl']);
            res.tplPathGen = checkForInput(varargin, 'tplPathGen', ...
                [res.templatePath, res.tplNameGen, '.tpl']);
            res.envFile = checkForInput(varargin, 'envFile', ...
                [res.dataFolder, filesep, res.envFilename]);
        end
        
        function initDir(obj)            
            if ~exist(obj.dataFolder, 'dir')
                mkdir(obj.dataFolder);
                disp(['Created ', obj.dataFolder, '.']);
            else
                disp([obj.dataFolder, ' already exists.']);
            end
            
            tpl = [obj.dataFolder, filesep, obj.tplName, '.tpl'];
            tplGen = [obj.dataFolder, filesep, obj.tplNameGen, '.tpl'];
            if ~exist(tpl, 'file')
                if exist(obj.tplPath, 'file')
                    copyfile(obj.tplPath, tpl);
                    disp(['Copied template ' obj.tplPath, ' to ', tpl, '.']);
                else
                    disp(['Template ', obj.tplPath, ' cannot be found.' ...
                        ' Please copy manually.']);
                end
            else
                disp(['Template file ', tpl, ' already exists.']);
            end
            if ~exist(tplGen, 'file')
                if exist(obj.tplPathGen, 'file')
                    copyfile(obj.tplPathGen, tplGen);
                    disp(['Copied template ' obj.tplPathGen, ' to ', ...
                        tplGen, '.']);
                else
                    disp(['Template ', obj.tplPathGen, ' cannot be found.' ...
                        ' Please copy manually.']);
                end
            else
                disp(['Template file ', tplGen, ' already exists.']);
            end
        end
        
        function makeEnv(obj, varargin)
            % Create a .mat file that contains the necessary environment
            % data. The varargin parameters that are checked depend on the
            % type of environments being created.
            % all:
            %     newEnv         Force overwriting old environment file.
            % envType = 'charge':
            %     mag            Maximum charge magnitude.
            %     cubSize        Cube size.
            %     cent           1 x 3 displacement vector for center of
            %                    molecule.
            % envType = 'field:'
            %     fieldType      Available field directions.
            %     fieldMethod    'random' or 'systematic'.
            %     magMin         Minimum field strength to use.
            %     magMax         Maximum field strength to use.
            % fieldMethod = 'systematic'
            %     magStep        Amount to increment field strength each
            %                    time.
            
            newEnv = checkForInput(varargin, 'newEnv', 0);
            
            % Avoid accidental overwrites.
            if (~newEnv && exist(obj.envFile, 'file'))
                disp(['Environments already exist. Pass in "''', ...
                    'newEnv''', ', 1" to overwrite.']);
                load(obj.envFile);
            else
                obj.nenv = checkForInput(varargin, 'nenv', 100);
                
                % Do charged environments.
                if strncmpi(obj.envType, 'charge', 6)
                    obj.makeChargeEnv(varargin);
                
                % Do fields.
                elseif strncmpi(obj.envType, 'field', 5)
                    obj.makeFieldEnv(varargin);
                    
                % Catch bad input.
                else
                    error('Invalid environment type.');
                end
            end
        end
        
        function makeChargeEnv(obj, varargin)
            mag = checkForInput(varargin, 'mag', 15.0);
            cubSize = checkForInput(varargin, 'cubSize', [6,6,6]);
            cent = checkForInput(varargin, 'cent', [0.77; 0; 0]);
            env = cell(1, obj.nenv);
            for ienv = 1:obj.nenv
                temp = Environment.newCube(cubSize,mag);
                temp.displace(cent);
                env{ienv} = temp;
            end
            save(obj.envfile,'env');
        end
        
        function makeFieldEnv(obj, varargin)
            fieldOpt = [ 1 0 0; 0 1 0; 0 0 1; 2 0 0; 0 2 0; 0 0 2; ...
                1 1 0; 1 0 1; 0 1 1; 3 0 0; 0 3 0; 0 0 3; 2 1 0; ...
                1 2 0; 2 0 1; 1 0 2; 0 2 1; 0 1 2; 1 1 1; 4 0 0; ...
                0 4 0; 0 0 4; 3 1 0; 1 3 0; 3 0 1; 1 0 3; 0 3 1; ...
                0 1 3; 2 2 0; 2 0 2; 0 2 2 ];
            fieldType = checkForInput(varargin, 'fieldType', fieldOpt);
            fieldMethod = checkForInput(varargin, 'fieldMethod', ...
                'random');
            magMin = checkForInput(varargin, 'magMin', 50);
            magMax = checkForInput(varargin, 'magMax', 800);
            
            % Randomly generated fields.
            % Only one field at a time is used right now.
            if strncmpi(fieldMethod, 'random', 6)
                env = cell(1, obj.nenv);
                for ifield = 1:obj.nenv
                    new = Environment;
                    new.nfield = 1;
                    tmp = 0;
                    while tmp == 0
                        tmp = int16(rand * size(fieldType,1));
                    end
                    new.fieldType = fieldType(tmp, :);
                    new.fieldMag = int16(rand * (magMin - magMax)) ...
                        + magMin;
                    new.ncharge = 0; new.rho = 0; new.r = 0;
                    env{ifield} = new;
                end
                
                % Systematically generated fields.
                % Only one field at a time is used right now.
            elseif strncmpi(fieldMethod, 'systematic', 10);
                magStep = checkForInput(varargin, 'magStep', 50);
                i = 1;
                for ifield = 1:size(fieldType, 1)
                    for imag = magMin:magStep:magMax
                        new = Environment;
                        new.nfield = 1;
                        new.fieldType = fieldType( ifield, : );
                        new.fieldMag = imag;
                        new.ncharge = 0; new.rho = 0; new.r = 0;
                        env{i} = new;
                        i = i + 1;
                    end
                end
            end
            save(obj.envFile, 'env');
        end
    
        function [HL, LL] = runData(obj, varargin)
            % Check if data already exists.
            if (exist(obj.datFile,'file'))
                disp('loading existing data');
                load(obj.datFile);
            else            
                % Load the environments.
                if (exist(obj.envFile, 'file'))
                    disp('loading existing environments');
                    load(obj.envFile, 'env');
                else
                    error('No environments found. Aborting...');
                end
                
                npar = size(obj.pars, 2);
                HL = cell(npar, length(obj.HLbasis));
                LL = cell(npar, 3);
                
                for ipar = 1:size(obj.pars,2)
                    par = obj.pars{ipar};
                    disp('pars: ');
                    disp(par);
                    
                    config = Fragment.defaultConfig();
                    config.method = obj.method;
                    config.par = par;
                    
                    % HL
                    for ihl = 1:size(obj.HLbasis, 2)
                        config.template = obj.tplName;
                        config.basisSet = obj.HLbasis{ihl};
                        disp(['ipar ',num2str(ipar),' loading HL ',num2str(ihl)]);
                        frag1 = Fragment(obj.dataRoot, config);
                        for ienv = 1:obj.nenv
                            display(['HL env ',num2str(ienv)]);
                            frag1.addEnv(env{ienv});
                        end
                        HL{ipar,ihl} = frag1;
                    end
                    % LL 1
                    config.basisSet = 'STO-3G';
                    frag2 = Fragment(obj.dataRoot, config);
                    disp(['ipar ',num2str(ipar),' loading LL 1']);
                    for ienv = 1:obj.nenv
                        display(['LL env ',num2str(ienv)]);
                        frag2.addEnv(env{ienv});
                    end
                    LL{ipar,1} = frag2;
                    
                    % LL 2
                    config.template = obj.tplNameGen;
                    config.basisSet = 'GEN';
                    config.par = [par 0.9 0.9 0.9 0.9 0.9];
                    frag3 = Fragment(obj.dataRoot, config);
                    disp(['ipar ',num2str(ipar),' loading LL 2']);
                    for ienv = 1:obj.nenv
                        display(['LL env ',num2str(ienv)]);
                        frag3.addEnv(env{ienv});
                    end
                    LL{ipar,2} = frag3;
                    % LL 3
                    config.template = obj.tplNameGen;
                    config.basisSet = 'GEN';
                    config.par = [par 1.05 1.05 1.05 1.05 1.05];
                    disp(['ipar ',num2str(ipar),' loading LL 3']);
                    frag4 = Fragment(obj.dataRoot, config);
                    for ienv = 1:obj.nenv
                        display(['LL env ',num2str(ienv)]);
                        frag4.addEnv(env{ienv});
                    end
                    LL{ipar,3} = frag4;
                end
                
                % since even loading all the files will take time, we'll dave everything
                save(obj.datFile);
            end
        end    
    end
end


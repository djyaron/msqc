classdef Fragment < handle
   %FRAGMENT Holds all data for a fragment and environments.
   %   After a Fragment is initialized, all necessary Gaussian calculations
   %   are performed and the results are stored.
   
   properties (SetAccess = private)
      config;        % structure holding calculation inputs
      dataPath;      % location of template and data files
      templateText;  % text from the template file
      gaussianFile;  % gaussian job (i.e. input) file (with charge keyword)
      fileprefix;    % prefix for files, without environment numbers
      
      natom   % number of atoms in fragment
      nelec   % number of electrons in the fragment
      Z       % (1,natom) atomic numbers of the atoms
      rcart   % (3,natom) cartesian coordinates of the atoms
      npar    % number of parameters in template file
      
      nbasis  % number of atomic (and molecular) basis functions
      H1      % (nbasis,nbasis) full H1 operator of fragment
      H1en;   % (nbasis,nbasis,natom) electron-nuclear interaction
      KE;     % (nbasis,nbasis) kinetic energy
      H2;     % (nbasis,nbasis,nbasis,nbasis) 2-elec interactions
      S;      % (nbasis,nbasis) overlap
      Hnuc;   % nuclear-nuclear interaction energy
      
      Ehf      % Hartree Fock energy
      MP2      % MP2 Energy
      %CorrE   % Correlation Energy (MP2-Ehf)
      Eorb     % (nbasis,1)      molecular orbital energies
      orb      % (nbasis,nbasis) molecular orbital coefficients
      dipole   % (3,1)   dipole moment of molecule
      mulliken % (1,natom)  mulliken charge on the atoms
      
      nenv    %  numer of environments
      env     % (1,nenv)             environments
      H1Env   % (nbasis,nbasis,nenv) H1 due to environment
      %  full H1 in environment = H1 + H1Env(:,:,ienv)
      HnucEnv % (1,nenv)             Hnuc in environment
      
      EhfEnv    % (1,nenv)        Hartree-Fock energy in env
      MP2Env    % (1,nenv)        MP2 energy in env
      %CorrEenv % (1,nenv)        Correlation energy in env (MP2env-HFenv)
      EorbEnv;  % (nbasis,nenv)   molecular orbital energies in env
      orbEnv;   % (nbasis,nbasis,nenv) molecular orbitals in env
      dipoleEnv % (3,nenv) dipole moment in the environment
      
      basisAtom    % (nbasis,1) atom # on which the function is centered
      basisType    % (nbasis,1) l quantum number: 0=s 1=p 2=d 3=d etc
      basisSubType % (nbasis,1) m quantum number: s=1 p=3 d=6 (cartesian)
      basisNprims  % number of primitives in this function
      basisPrims   % {nbasis,1} cell array of matrices of size (2,nprims)
      %    with (1,:) being contraction coefficients and
      %         (2,:) being primimitive exponents
      
   end
   properties
      % TODO Need to do something about these
      gaussianPath = 'c:\g09w';
      gaussianExe  = 'g09.exe';
   end
   methods (Access = private)
      initializeData(obj);
   end
   methods (Static)
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
         res.par      = [];
      end
      [found,fileprefix] = findCalc(dataPath,config)
      [CorrE, MP2, Ehf, Eorb, orb, Nelectrons,  Z, rcart, ...
         dipole, mulliken, ...
         atom, type, subtype, nprims, prims ] = readfchk(fid1)
      [Eorb, orb, atom, Nelectrons, Ehf] = oldreadfchk(fid1)
      [S, H1, KE, H2, Enuc] = readpolyatom(fid1)
      [HL, LL] = dataMerge(datFiles, envs, saveFilename, varargin)
   end
   methods
      function res = Fragment(dataPathIn, configIn, useCache)
         %  dataPath = directory (including c:\ etc) for data storage
         %               do not include a \ at end of paths
         %               [defaults to 'data']
         %  configIn = configuration structure
         %               [defaults to 'Fragment.defaultConfig();
         %  useCache = look in dataPathIn for cached result, and save 
         %             new results in dataPathIn [defaults to true]
         if (nargin < 1)
            res.dataPath = 'data';
         else
            res.dataPath = dataPathIn;
         end
         if (nargin < 2)
            res.config = Fragment.defaultConfig();
         else
            res.config = configIn;
         end
         if (nargin < 3)
            useCache = true;
         end
         newline = char(10);
         if (useCache)
            [found,res.fileprefix] = ...
               Fragment.findCalc(res.dataPath,res.config);
         else
            found = false;
            res.fileprefix = '';
         end
         % Backword compatibility (load mat file)
         if (found && exist([res.fileprefix,'_calc.mat'],'file'))
            ftemp = [res.fileprefix,'_calc.mat'];
            prefixsave = res.fileprefix;
            dataPathsave = res.dataPath;
            load(ftemp, 'resFile' );
            res = resFile;
            res.fileprefix = prefixsave;
            res.dataPath = dataPathsave;
         else % load zip file, with generation if needed
            % If the cfg file exists, but not the zip file, delete the cfg
            % and proceed as though it was not found.
            if (found && ~exist([res.fileprefix, '.zip'], 'file'))
               delete([res.fileprefix, '_cfg.mat']);
               found = false;
            end
            
            res.templateText = fileread([res.dataPath,filesep,...
               res.config.template,'.tpl']);
            res.natom = size( strfind(res.templateText, 'ATOM'), 2);
            res.npar = size( strfind(res.templateText, 'PAR'), 2);
            
            basisSet = res.config.basisSet;
            method   = res.config.method;
            charge   = res.config.charge;
            spin     = res.config.spin;
            par      = res.config.par;
            
            % ___________________________________________________________
            % header for the Gaussian job file (input file)
            % Note: for single atom calcs below, 'scf=conventional' is replaced
            %       so if this keyword in header is changed, it needs to be 
            %       changed there as well
            header = ['%rwf=temp.rwf',newline,...
               '%nosave',newline,...
               '%chk=temp.chk',newline,...
               '# ',method,'/',basisSet, newline...
               'nosymm int=noraff iop(99/6=1) ',...
               'scf=conventional',' symm=noint', newline, newline, ...
               'title', newline,newline];
            
            % ____________________________________________________________
            % Create Scratch directory within g09 scratch directory, to do work
            
            % ctext will hold the Gaussian job file (input file)
            % begin with the above header text
            ctext = header;
            % charge and spin come next
            ctext = [ctext, num2str(charge), ' ', num2str(spin), newline];
            % For molecule specification, we first replace all ATOM# with spaces
            t1 = res.templateText;
            % Iterate in reverse order, or replacements will not work properly 
            % with more than 10 atoms.
            for iatom = res.natom:-1:1
               t1 = strrep(t1, ['ATOM',num2str(iatom)], ' ');
            end
            % And replace all PAR# with the parameter values.
            for ipar = res.npar:-1:1
               t1 = strrep(t1, ['PAR',num2str(ipar)], num2str(par(ipar),'%23.12f'));
            end
            ctext = [ctext, t1];
            
            res.gaussianFile = ctext;
            
            nparIn = size(res.config.par,1) * size(res.config.par,2);
            if (nparIn ~= res.npar)
               error(['template has ',num2str(res.npar),' parameters',...
                  ' while config contains ',num2str(nparIn),' pars']);
            end
            if (~found)
               temp1 = tempname('a'); % makes "a\uniquestring"
               uniqueStr = temp1(3:end);
               res.fileprefix = [res.dataPath,filesep,res.config.template, ...
                  '_',uniqueStr];
               % save config file, in *.mat format
               Cfile = res.config;
               save([res.fileprefix,'_cfg.mat'],  'Cfile' );
               % create and save zip file
               zipFile = [res.fileprefix,'.zip'];
               res.initializeZipData(zipFile);
            end
            zipFile = [res.fileprefix,'.zip'];
            res.loadZipData(zipFile);
             if (~useCache)
                delete(zipFile);
             end
         end
         res.nenv = 0;
         % Set the environment array to have the correct class type
         res.env = Environment.empty(0,0);
      end
      function setEnvSize(obj,nenvIn)
         clear obj.env;
         obj.env(1,nenvIn) = Environment;
         obj.H1Env = zeros(obj.nbasis, obj.nbasis, nenvIn);
         obj.EhfEnv = zeros(1,nenvIn);
         obj.MP2Env = zeros(1,nenvIn);
         obj.CorrEenv  = zeros(1,nenvIn);
         obj.EorbEnv = zeros(obj.nbasis, nenvIn);
         obj.orbEnv  = zeros(obj.nbasis,obj.nbasis,nenvIn);
         obj.dipoleEnv = zeros(3,nenvIn);
      end
   end % methods
end %

classdef Fragment < handle
   properties (SetAccess = private)      
      config;        % structure holding calculation inputs
      dataPath;      % location of template and data files
      templateText;  % text from the template file
      gaussianFile;  % gaussian job (i.e. input) file (with charge keyword)
      fileprefix;    % prefix for files, without environment numbers
      
      natom   % number of atoms in fragment
      npar    % number of parameters in template file
      nbasis  % number of atomic basis functions
      nelec   % number of electrons in the fragment
      H1      % (nbasis,nbasis) full H1 operator of fragment
      H1en;   % (nbasis,nbasis,natom) electron-nuclear interaction
      KE;   % (nbasis,nbasis) kinetic energy
      H2;     % (nbasis,nbasis,nbasis,nbasis) 2-elec interactions
      S;      % (nbasis,nbasis) overlap
      Hnuc;   % nuclear-nuclear interaction energy
      basisAtom % (nbasis,1) basisAtom(ibasis) = atom on which ibasis resides
      
      Ehf     % Hartree Fock energy
      Eorb    % (nbasis,1)      molecular orbital energies
      orb     % (nbasis,nbasis) molecular orbital coefficients
      
      nenv    %  numer of environments
      env     % (1,nenv)             environments
      H1Env   % (nbasis,nbasis,nenv) H1 due to environment
              %  full H1 in environment = H1 + H1Env(:,:,ienv)
      HnucEnv % (1,nenv)             Hnuc in environment

      EhfEnv   % (1,nenv)        Hartree-Fock energy in env
      EorbEnv; % (nbasis,nenv)   molecular orbital energies in env
      orbEnv;  % (nbasis,nbasis,nenv) molecular orbitals in env  
   end
   properties
      % TODO Need to do something about these
      gaussianPath = 'c:\g03w';
      gaussianExe  = 'g03.exe';
   end
   methods (Access = private)
      initializeData(obj);
      initializeDataInEnv(obj,env);
   end
   methods (Static)
      function res = defaultConfig()
         %  templateFile = name of template file, in dataPathIn
         %               [defaults to 'template.txt']
         %  basisSet = basis set keyword (Gaussian format)
         %               [defaults to 'STO-3G']
         %  charge = charge on the fragment
         %             [defaults to 0]
         %  spin   = spin (multiplicity) of the fragment,
         %             using Gaussian convention
         %             [defaults to 1]
         res.template = 'template';
         res.basisSet = 'STO-3G';
         res.charge   = 0;
         res.spin     = 1;
         res.par      = [];
      end
      [found,fileprefix] = findCalc(dataPath,config)
      [Eorb, orb, atom, Nelectrons, Ehf] = readfchk(fid1)
      [S, H1, KE, H2, Enuc] = readpolyatom(fid1)
   end
   methods
      function res = Fragment(dataPathIn, configIn)
         %  dataPath = directory (including c:\ etc) for data storage
         %               do not include a \ at end of paths
         %               [defaults to 'data']
         %  configIn = configuration structure
         %               [defaults to 'Fragment.defaultConfig();
         if (nargin < 1)
            res.dataPath = 'data';
         else
            res.dataPath = dataPathIn;
         end
         if (nargin < 2)
            res.config = Fragment.defaultConfig();
         else
            res.config = configIn;
            [found,res.fileprefix] = ...
               Fragment.findCalc(res.dataPath,res.config);
            if (found)
               ftemp = [res.fileprefix,'_calc.mat'];
               load(ftemp, 'resFile' );
               res = resFile;
            else
               res.templateText = fileread([res.dataPath,'\',...
                  res.config.template,'.tpl']);
               res.natom = size( strfind(res.templateText, 'ATOM'), 2);
               res.npar = size( strfind(res.templateText, 'PAR'), 2);
               res.initializeData();
               resFile = res;
               Cfile = res.config;
               save([res.fileprefix,'_cfg.mat'],  'Cfile' );
               save([res.fileprefix,'_calc.mat'], 'resFile' );
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
         obj.EorbEnv = zeros(obj.nbasis, nenvIn);
         obj.orbEnv  = zeros(obj.nbasis,obj.nbasis,nenvIn);
      end
   end % methods
end %
classdef Model1 < handle
   properties (SetAccess = private)            
      natom   % number of atoms in fragment
      nbasis  % number of atomic basis functions
      nelec   % number of electrons in the fragment
      nenv    %  numer of environments
      env     % (1,nenv)             environments
      
      Ehf     % Hartree Fock energy
      Eorb    % (nbasis,1)      molecular orbital energies
      orb     % (nbasis,nbasis) molecular orbital coefficients
      
      EhfEnv   % (1,nenv)        Hartree-Fock energy in env
      EorbEnv; % (nbasis,nenv)   molecular orbital energies in env
      orbEnv;  % (nbasis,nbasis,nenv) molecular orbitals in env  

   end      
   properties
      H1      % (nbasis,nbasis) full H1 operator of fragment
      H1en;   % (nbasis,nbasis,natom) electron-nuclear interaction
      KE;     % (nbasis,nbasis) kinetic energy
      H2;     % (nbasis,nbasis,nbasis,nbasis) 2-elec interactions
      S;      % (nbasis,nbasis) overlap
      Hnuc;   % nuclear-nuclear interaction energy
      basisAtom % (nbasis,1) basisAtom(ibasis) = atom on which ibasis resides
      
      H1Env   % (nbasis,nbasis,nenv) H1 in environments
      HnucEnv % (1,nenv)             Hnuc in environment

   end
   methods (Static)
      [orb,Eorb,Ehf] = hartreeFock(frag,env,epsIn)
   end
   methods
      function res = Model1(frag)
         res.natom      = frag.natom;
         res.nbasis     = frag.nbasis;
         res.nelec      = frag.nelec;
         res.nenv       = frag.nenv;
         res.env        = frag.env;
         res.H1         = frag.H1;
         res.H1en       = frag.H1en;
         res.KE         = frag.KE;
         res.H2         = frag.H2;
         res.S          = frag.S;
         res.Hnuc       = frag.Hnuc;
         res.basisAtom  = frag.basisAtom;
         res.H1Env      = frag.H1Env;
         res.HnucEnv    = frag.HnucEnv;
      end
   end % methods
end %
classdef Model2 < handle
%{
   Interpolate between narrow and diffuse STO-3G matrix elements. 
   Hartree-Fock routine needs the following methods:
      H1 H1env Hnuc HnucEnv H2 S nelec nbasis
%}
   properties
      frag;  % Fragment with regular STO-3G basis set
      fnar;  % Fragment with narrow  STO-3G basis set
      fdif;  % Fragment with diffuse STO-3G basis set
      fHL;   % High level data
      
      Ehf     % Hartree Fock energy
      Eorb    % (nbasis,1)      molecular orbital energies
      orb     % (nbasis,nbasis) molecular orbital coefficients
      EhfEnv   % (1,nenv)        Hartree-Fock energy in env
      EorbEnv; % (nbasis,nenv)   molecular orbital energies in env
      orbEnv;  % (nbasis,nbasis,nenv) molecular orbitals in env
      
      
      natom   % number of atoms in fragment
      nelec   % number of electrons in the fragment
      Z       % (1,natom) atomic numbers of the molecules
      rcart   % (3,natom) cartesian coordinates of the atoms
      nenv    

      nbasis  % number of atomic (and molecular) basis functions
      basisAtom  % (nbasis,1) atom # on which the function is centered 
      basisType  % (nbasis,1) l quantum number: 0=s 1=p 2=d 3=d etc
      basisSubType % (nbasis,1) m quantum number: s=1 p=3 d=6 (cartesian)
      onAtom     % {natom,1}  list of basis functions on iatom
      
      par     % current list of parameters for fitting routine
      E2tmp   % temporary store of E2 in environments, for errApprox
   end
   methods (Static)
      function res = mix(x, v1, v2)
         % mix objects v1 and v2, using parameter x. 
         %   for x << 0, we get v1, and x>>0 we get v2, with the 
         %   switch from v1 to v2 occuring mostly as x=-1..1
         c1 = (tanh(x)+1)/2.0;
         c2 = 1-c1;
         res = c2 * v1 + c1 * v2;
      end
   end
   methods
      function res = Model2(frag_,fnar_, fdif_,fHL_)
         res.frag = frag_;
         res.fnar = fnar_;
         res.fdif = fdif_;
         res.fHL  = fHL_;
         res.natom = frag_.natom;
         res.nelec = frag_.nelec;
         res.Z     = frag_.Z;
         res.rcart = frag_.rcart;
         res.nenv  = frag_.nenv;
         res.nbasis = frag_.nbasis;
         res.basisAtom = frag_.basisAtom;
         res.basisType = frag_.basisType;
         res.basisSubType = frag_.basisSubType;
         for iatom = 1:res.natom
            res.onAtom{iatom,1} = find(res.basisAtom == iatom);
         end
      end
      function res = updateDensity(obj,par)
          obj.par = par;
          obj.solveHF;
          obj.E2tmp = zeros(1,obj.nenv);
          for ienv = 1:obj.nenv
              obj.E2tmp(ienv) = sum(sum(sum(sum(obj.partitionE2(ienv)))));
          end
      end
      function res = errApprox(obj,par)
          obj.par = par;
          res = zeros(1,obj.nenv);
          for ienv=1:obj.nenv
              E1 = sum(sum(obj.partitionE1(ienv)));
              E2 = obj.E2tmp(ienv);
              Enuc = obj.Hnuc(ienv);
              res(ienv) = (E1+E2+Enuc) - obj.fHL.EhfEnv(ienv);
          end
      end
      
      function res = err(obj, par)
         obj.par = par;
         obj.solveHF;
         ic = 0;
         res = zeros(1,obj.nenv);
         %res(ic) = obj.Ehf - obj.fHL.Ehf;
         for i=1:obj.nenv
             ehl = obj.fHL.EhfEnv(i);
             emod = obj.EhfEnv(i);
             ic = ic+1;
             res(ic) = emod - ehl;
         end
         %den = obj.density;
         %denHL = obj.fHL.density;
      end
      function res = H1(obj, ienv)
         % parameters for the H1 diagonal energies
         diagPar = zeros(6,1);
         par = obj.par;
         diagPar(1,1) = par(1); % diagonal for H
         diagPar(6,1) = par(2); % diagonal for C
         % Bonding parameters between atoms
         bondPar = zeros(6,6);
         bondPar(6,6) = par(3); % C C bonds
         bondPar(1,6) = par(4); % C H bonds
         bondPar(6,1) = bondPar(1,6);
         res   = obj.frag.H1;
         for iatom = 1:obj.natom
            ran = obj.onAtom{iatom}; % range of functions on iatom
            % substract off current value
            res(ran,ran) = res(ran,ran) - obj.frag.H1en(ran,ran,iatom)- obj.frag.KE(ran,ran,iatom);
            % add on scaled value
            res(ran,ran) = res(ran,ran) + Model2.mix( diagPar( obj.Z(iatom) ), ...
                obj.fnar(ran,ran,iatom), obj.fdif(ran,ran,iatom));
         end
         bond = bonded();
         for iatom = 1:obj.natom
             for jatom = 1:obj.natom
%              distAtom(iatom,jatom) = sqrt(((cart(1,iatom)-(cart(1,jatom)))^2)...
%              +((cart(2,iatom)-(cart(2,jatom)))^2)+((cart(3,iatom)-(cart(3,jatom)))^2));
             if  bond(iatom,jatom) ==1
                 iran = obj.onAtom{iatom};
                 jran = obj.onAtom{jatom};
                 res(iran,jran) = res(iran,jran) ... 
                  - obj.frag.H1en(iran,jran,iatom) ...
                  - obj.frag.H1en(iran,jran,jatom) ...
                  - obj.frag.KE(iran,jran);
               % add on scaled value
               vnar = obj.fnar.H1en(iran,jran,iatom) ...
                  + obj.fnar.H1en(iran,jran,jatom) ...
                  + obj.fnar.KE(iran,jran);
               vdif = obj.fdif.H1en(iran,jran,iatom) ...
                  + obj.fdif.H1en(iran,jran,jatom) ...
                  + obj.fnar.KE(iran,jran);
               res(iran,jran) = res(iran,jran) + ...
                  obj.mix( bondPar(obj.Z(iatom),obj.Z(jatom)), vnar, vdif);
             end
             end
         end
         if (ienv > 0)
            res = res + obj.frag.H1Env(:,:,ienv);
         end
      end
      function res = Hnuc(obj,ienv)
         if (ienv == 0)
            res = obj.frag.Hnuc;
         else
            res = obj.frag.HnucEnv(ienv);
         end
      end
      function res = H2(obj)
         res = obj.frag.H2;
      end
      function res = S(obj)
         res = obj.frag.S;
      end
   end % methods
end %
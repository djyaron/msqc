classdef Model2 < handle
%{
   Interpolate between narrow and diffuse STO-3G matrix elements. 
%}
   properties
      % Input to the class
      frag;  % Fragment with regular STO-3G basis set
      fnar;  % Fragment with narrow  STO-3G basis set
      fdif;  % Fragment with diffuse STO-3G basis set
      
      % Most recent predictions of the model
      Ehf     % Hartree Fock energy
      Eorb    % (nbasis,1)      molecular orbital energies
      orb     % (nbasis,nbasis) molecular orbital coefficients
      EhfEnv   % (1,nenv)        Hartree-Fock energy in env
      EorbEnv; % (nbasis,nenv)   molecular orbital energies in env
      orbEnv;  % (nbasis,nbasis,nenv) molecular orbitals in env
      
      % Useful properties initialized in the constructor
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
      valAtom    % {natom,type} list of valence basis functions of type 
                 %              (1-s 2-p) on iatom
      isBonded   % (natom,natom)  1 if atoms are bonded, 0 otherwise
      
      % Current parameters, and some temporary values
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
      function res = Model2(frag_,fnar_, fdif_)
         res.frag = frag_;
         res.fnar = fnar_;
         res.fdif = fdif_;
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
         for iatom = 1:res.natom
            % kind of a hack. For s orbitals, we take the maximum value
            % from the list of basis functions on the atom, since this 
            % will be the valence orbital (2s instead of 1s for C)
            res.valAtom{iatom,1} = find(res.basisAtom == iatom & ...
                                            res.basisType == 0, 1, 'last' );
            % For p orbitals, we just take the ones that matach
            res.valAtom{iatom,2} = find(res.basisAtom == iatom & ...
                                              res.basisType == 1);
         end
         res.isBonded = zeros(res.natom,res.natom);
         for iatom = 1:res.natom
            for jatom = 1:res.natom
               res.isBonded(iatom,jatom) = res.bonded(iatom,jatom);
            end
         end
      end
      function res = updateDensity(obj,par)
         % updates density for use in obj.errApprox, including saving the
         % two electron energy. 
          obj.par = par;
          obj.solveHF;
          obj.E2tmp = zeros(1,obj.nenv);
          for ienv = 1:obj.nenv
              obj.E2tmp(ienv) = sum(sum(sum(sum(obj.partitionE2(ienv)))));
          end
      end
      function res = H1check(obj, ienv)
         % parameters for the H1 diagonal energies
         diagParKE = zeros(6,2); % element and type(s,p)
         diagparEN = zeros(6,2);
         par = obj.par;
         diagParKE(1,1) = par(1); % diagonal for H
         diagParKE(6,1) = par(2); % s diagonal for C
         diagParKE(6,2) = par(3); % s diagonal for C
         diagParEN(1,1) = par(4); % diagonal for H
         diagParEN(6,1) = par(5); % s diagonal for C
         diagParEN(6,2) = par(6); % s diagonal for C
         % Bonding parameters between atoms
         bondParKE = zeros(6,2,6,2);
         bondParKE(1,1,6,1) = par(7); % H Cs bonds
         bondParKE(6,1,1,1) = par(7); %
         bondParKE(1,1,6,2) = par(8); % H Cp bonds
         bondParKE(6,2,1,1) = par(8); %
         bondParKE(6,1,6,1) = par(9); % Cs Cs bonds
         bondParKE(6,1,6,2) = par(10); % Cs Cp bonds
         bondParKE(6,2,6,1) = par(10); %
         bondParKE(6,2,6,2) = par(11); % Cp Cp bonds
         bondParEN = zeros(6,2,6,2);
         bondParEN(1,1,6,1) = par(12); % H Cs bonds
         bondParEN(6,1,1,1) = par(12); %
         bondParEN(1,1,6,2) = par(13); % H Cp bonds
         bondParEN(6,2,1,1) = par(13); %
         bondParEN(6,1,6,1) = par(14); % Cs Cs bonds
         bondParEN(6,1,6,2) = par(15); % Cs Cp bonds
         bondParEN(6,2,6,1) = par(15); %
         bondParEN(6,2,6,2) = par(16); % Cp Cp bonds
         
         % start with H1 matrix of unmodified STO-3G
         res   = obj.frag.H1;
         
         % On atom modifications: KE and H1en with atom on which orbitals
         %  reside. 
         for iatom = 1:obj.natom
            if (obj.Z(iatom) == 1) % Hydrogen
               maxType = 1;
            else % heavy atom
               maxType = 2;
            end
            for itype = 1:maxType
               ran = obj.valAtom{iatom,itype}; % range of valence orbs (1s)
               res(ran,ran) = res(ran,ran) - obj.frag.KE(ran,ran) ...
                  + Model2.mix( diagParKE(obj.frag.Z(iatom),itype), ...
                  obj.fnar.KE(ran,ran), ...
                  obj.fdif.KE(ran,ran) );
               res(ran,ran) = res(ran,ran) - obj.frag.H1en(ran,ran,iatom) ...
                  + Model2.mix( diagParEN(obj.frag.Z(iatom),itype), ...
                  obj.fnar.H1en(ran,ran,iatom),...
                  obj.fdif.H1en(ran,ran,iatom));
            end
         end
         % Bonding terms 
         for iatom = 1:obj.natom
            if (obj.Z(iatom) == 1) % Hydrogen
               imaxType = 1;
            else % heavy atom
               imaxType = 2;
            end
            for jatom = 1:obj.natom
               if (obj.Z(jatom) == 1) % Hydrogen
                  jmaxType = 1;
               else % heavy atom
                  jmaxType = 2;
               end
               if (obj.isBonded(iatom,jatom))
                  for itype = 1:imaxType
                     for jtype = 1:jmaxType
                        iran = obj.valAtom{iatom,itype};
                        jran = obj.valAtom{jatom,jtype};
                        iZ = obj.frag.Z(iatom);
                        jZ = obj.frag.Z(jatom);
                        % modify kinetic energy
                        res(iran,jran) = res(iran,jran) ...
                           - obj.frag.KE(iran,jran) ...
                           + obj.mix(bondParKE(iZ,itype,jZ,jtype), ...
                           obj.fnar.KE(iran,jran), ...
                           obj.fdif.KE(iran,jran));
                        % modify elec-nuc interaction
                        res(iran,jran) = res(iran,jran) ...
                           - obj.frag.H1en(iran,jran,iatom) ...
                           - obj.frag.H1en(iran,jran,jatom);
                        % add on scaled value
                        vnar = obj.fnar.H1en(iran,jran,iatom) ...
                           + obj.fnar.H1en(iran,jran,jatom);
                        vdif = obj.fdif.H1en(iran,jran,iatom) ...
                           + obj.fdif.H1en(iran,jran,jatom);
                        res(iran,jran) = res(iran,jran) + ...
                           obj.mix(bondParEN(iZ,itype,jZ,jtype),...
                           vnar, vdif);
                     end
                  end
               end
            end
         end
         if (ienv > 0)
            res = res + obj.frag.H1Env(:,:,ienv);
         end
      end
      function res = H1(obj, ienv)
         if (nargin < 2)
            ienv = 0;
         end
         res = obj.KE;
         for iatom = 1:obj.natom
            res = res + obj.H1en(iatom);
         end
         if (ienv > 0)
            res = res + obj.frag.H1Env(:,:,ienv);
         end
      end
      function res = KE(obj)
         par = obj.par;
         diagParKE = zeros(6,2); % element and type(s,p)
         diagParKE(1,1) = par(1); % diagonal for H
         diagParKE(6,1) = par(2); % s diagonal for C
         diagParKE(6,2) = par(3); % s diagonal for C
         % Bonding parameters between atoms
         bondParKE = zeros(6,2,6,2);
         bondParKE(1,1,6,1) = par(7); % H Cs bonds
         bondParKE(6,1,1,1) = par(7); %
         bondParKE(1,1,6,2) = par(8); % H Cp bonds
         bondParKE(6,2,1,1) = par(8); %
         bondParKE(6,1,6,1) = par(9); % Cs Cs bonds
         bondParKE(6,1,6,2) = par(10); % Cs Cp bonds
         bondParKE(6,2,6,1) = par(10); %
         bondParKE(6,2,6,2) = par(11); % Cp Cp bonds
         
         % start with H1 matrix of unmodified STO-3G
         res   = obj.frag.KE;
         
         % On atom modifications: KE and H1en with atom on which orbitals
         %  reside. 
         for iatom = 1:obj.natom
            if (obj.Z(iatom) == 1) % Hydrogen
               maxType = 1;
            else % heavy atom
               maxType = 2;
            end
            for itype = 1:maxType
               ran = obj.valAtom{iatom,itype}; % range of valence orbs (1s)
               res(ran,ran) = res(ran,ran) - obj.frag.KE(ran,ran) ...
                  + Model2.mix( diagParKE(obj.frag.Z(iatom),itype), ...
                  obj.fnar.KE(ran,ran), ...
                  obj.fdif.KE(ran,ran) );
            end
         end
         % Bonding terms 
         for iatom = 1:obj.natom
            if (obj.Z(iatom) == 1) % Hydrogen
               imaxType = 1;
            else % heavy atom
               imaxType = 2;
            end
            for jatom = 1:obj.natom
               if (obj.Z(jatom) == 1) % Hydrogen
                  jmaxType = 1;
               else % heavy atom
                  jmaxType = 2;
               end
               if (obj.isBonded(iatom,jatom))
                  for itype = 1:imaxType
                     for jtype = 1:jmaxType
                        iran = obj.valAtom{iatom,itype};
                        jran = obj.valAtom{jatom,jtype};
                        iZ = obj.frag.Z(iatom);
                        jZ = obj.frag.Z(jatom);
                        % modify kinetic energy
                        res(iran,jran) = res(iran,jran) ...
                           - obj.frag.KE(iran,jran) ...
                           + obj.mix(bondParKE(iZ,itype,jZ,jtype), ...
                           obj.fnar.KE(iran,jran), ...
                           obj.fdif.KE(iran,jran));
                     end
                  end
               end
            end
         end
      end
      function res = H1en(obj, iatom)
         % parameters for the H1 diagonal energies
         diagparEN = zeros(6,2);
         par = obj.par;
         diagParEN(1,1) = par(4); % diagonal for H
         diagParEN(6,1) = par(5); % s diagonal for C
         diagParEN(6,2) = par(6); % s diagonal for C
         % Bonding parameters between atoms
         bondParEN = zeros(6,2,6,2);
         bondParEN(1,1,6,1) = par(12); % H Cs bonds
         bondParEN(6,1,1,1) = par(12); %
         bondParEN(1,1,6,2) = par(13); % H Cp bonds
         bondParEN(6,2,1,1) = par(13); %
         bondParEN(6,1,6,1) = par(14); % Cs Cs bonds
         bondParEN(6,1,6,2) = par(15); % Cs Cp bonds
         bondParEN(6,2,6,1) = par(15); %
         bondParEN(6,2,6,2) = par(16); % Cp Cp bonds
         
         % start with H1 matrix of unmodified STO-3G
         res   = obj.frag.H1en(:,:,iatom);
         if (obj.Z(iatom) == 1) % Hydrogen
            imaxType = 1;
         else % heavy atom
            imaxType = 2;
         end
         for itype = 1:imaxType
            ran = obj.valAtom{iatom,itype}; % range of valence orbs (1s)
            res(ran,ran) = res(ran,ran) - obj.frag.H1en(ran,ran,iatom) ...
               + Model2.mix( diagParEN(obj.frag.Z(iatom),itype), ...
               obj.fnar.H1en(ran,ran,iatom),...
               obj.fdif.H1en(ran,ran,iatom));
         end
         % Bonding terms 
         for jatom = 1:obj.natom
            if (obj.isBonded(iatom,jatom))
               if (obj.Z(jatom) == 1) % Hydrogen
                  jmaxType = 1;
               else % heavy atom
                  jmaxType = 2;
               end
               for itype = 1:imaxType
                  for jtype = 1:jmaxType
                     iran = obj.valAtom{iatom,itype};
                     jran = obj.valAtom{jatom,jtype};
                     iZ = obj.frag.Z(iatom);
                     jZ = obj.frag.Z(jatom);
                     % modify elec-nuc interaction
                     res(iran,jran) = res(iran,jran) ...
                        - obj.frag.H1en(iran,jran,iatom);
                     % add on scaled value
                     vnar = obj.fnar.H1en(iran,jran,iatom);
                     vdif = obj.fdif.H1en(iran,jran,iatom);
                     res(iran,jran) = res(iran,jran) + ...
                        obj.mix(bondParEN(iZ,itype,jZ,jtype),...
                        vnar, vdif);
                     
                     res(jran,iran) = res(jran,iran) ...
                        - obj.frag.H1en(jran,iran,iatom);
                     % add on scaled value
                     vnar = obj.fnar.H1en(jran,iran,iatom);
                     vdif = obj.fdif.H1en(jran,iran,iatom);
                     res(jran,iran) = res(jran,iran) + ...
                        obj.mix(bondParEN(iZ,itype,jZ,jtype),...
                        vnar, vdif);
                  end
               end
            end
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
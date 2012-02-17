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
      
      % model characteristics
      sepKE  % if 1, use different parameters for KE and H1en
      sepSP  % if 1, use different parameters for s and p orbs
      rhodep % if 1, parameters have a charge dependence
      mixType % if 0, use tanh; if 1, use linear
      
      % Current parameters
      par   % current list of parameters for fitting routine
            % KE diagonal    1: H  2: Cs  3: Cp
            % Hen diagonal   4: H  5: Cs  6: Cp
            % KE bonding     7: H-Cs  8: H-Cp  9: Cs-Cs 10: Cs-Cp 11: Cp-Cp  
            % Hen bonding   12: H-Cs 13: H-Cp 14: Cs-Cs 15: Cs-Cp 16: Cp-Cp
            % KE charge     17: H    18:Cs   19:Cp
      rhozero % (natom,nenv+1) mulliken charges of frag_ 
   end
   properties (Transient)
      densitySave   % cell array {1:nenv+1} of most recent density matrices 
                    % used to start HF iterations
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
         res.sepKE = 1;
         res.sepSP = 1;
         res.rhodep = 1;
         res.rhozero = zeros(frag_.natom,frag_.nenv+1);
         for ienv = 0:frag_.nenv
            res.rhozero(:,ienv+1) = frag_.mcharge(ienv)';
         end
         res.mixType = 0;
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
         res.densitySave = cell(1,res.nenv+1);
      end
      function res = mcharge(obj, ienv)
          % Calculates the mulliken charges on each atom
          Q = zeros(size(obj.Z));
          GAP = zeros(1,obj.natom);
          P = obj.density(ienv).*obj.S;
          GOP = sum(P,1);
          arange = cell(obj.natom,1);
          for iatom = 1:obj.natom
              arange{iatom} = find(obj.basisAtom == iatom);
          end
          for i = 1:obj.natom
              GAP(i) = sum(GOP(1,arange{i}));
              Q(i) = obj.Z(i)-GAP(i);
          end
          res = Q;
      end
      function res = mix(obj, x, v1, v2)
         if (obj.mixType == 0)
            % mix objects v1 and v2, using parameter x.
            %   for x << 0, we get v1, and x>>0 we get v2, with the
            %   switch from v1 to v2 occuring mostly as x=-1..1
            c1 = (tanh(x)+1)/2.0;
            c2 = 1-c1;
            res = c2 * v1 + c1 * v2;
         else
            % want linear mix, with (v1+v2)/2 when x=0
            % res = (v1+v2)/2 + x (v2-v1)/2
            % The bounds are: res = v1 at x = -1;
            %                 res = v2 at x = 1;
            % potentially faster (since v's are matrices while x is scalar)
            res = ((1.0-x)/2.0) * v1 + ((1.0+x)/2.0) * v2;
         end
      end
      function res = npar(obj)
         if (obj.sepKE && obj.sepSP && obj.rhodep)
            res = 19;
         end
         if (obj.sepKE && obj.sepSP && ~obj.rhodep)
            res = 16;
         end
         if (~obj.sepKE && ~obj.sepSP)
            res = 4;
         end
         if (~obj.sepKE && obj.sepSP)
            res = 8;
         end
         if (obj.sepKE && ~obj.sepSP)
            res = 8;
         end
      end    
      function res = setPar(obj,pIn)
         if (obj.sepKE && obj.sepSP)
            obj.par = pIn;
         end
         if (~obj.sepKE && ~obj.sepSP)
            p = zeros(1,16);
            p([1 4])              = pIn(1); % H diag
            p([2 3 5 6])          = pIn(2); % C diag
            p([7 8 12 13])        = pIn(3); % C H bond
            p([9 10 11 14 15 16]) = pIn(4); % C C bond
            obj.par = p;
         end
         if (~obj.sepKE && obj.sepSP)
            p = zeros(1,16);
            p([1 4])              = pIn(1); % H diag
            p([2 5])              = pIn(2); % Cs diag
            p([3 6])              = pIn(3); % Cp diag
            p([7 12])             = pIn(4); % H Cs bond
            p([8 13])             = pIn(5); % H Cp bond
            p([9 14])             = pIn(6); % Cs Cs bond
            p([10 15])            = pIn(7); % Cs Cp bond
            p([11 16])            = pIn(8); % Cp Cp bond
            obj.par = p;
         end
         if (obj.sepKE && ~obj.sepSP)
            p = zeros(1,16);
            p(1)                  = pIn(1); % H diag KE
            p(4)                  = pIn(2); % H diag EN
            p([2 3])              = pIn(3); % C diag KE
            p([5 6])              = pIn(4); % C diag EN
            p([7 8])              = pIn(5); % C H bond KE
            p([12 13])            = pIn(6); % C H bond EN
            p([9 10 11])          = pIn(7); % C C bond KE
            p([14 15 16])         = pIn(8); % C C bond EN
            obj.par = p;
         end
         res = obj.par;
      end               
      function res = mapPar(obj,pFull)
         % Given a full parameter vector, set guesses for the next level
         if (obj.sepKE && obj.sepSP)
            res = pFull;
         end
         if (~obj.sepKE && ~obj.sepSP)
            res = zeros(1,4);
            res(1) = pFull(1);
            res(2) = pFull(2);
            res(3) = pFull(7);
            res(4) = pFull(9);
         end
         if (~obj.sepKE && obj.sepSP)
            res = zeros(1,8);
            res(1) = pFull(1);
            res(2) = pFull(2);
            res(3) = pFull(3);
            res(4) = pFull(7);
            res(5) = pFull(8);
            res(6) = pFull(9);
            res(7) = pFull(10);
            res(8) = pFull(11);
         end
         if (obj.sepKE && ~obj.sepSP)
            res = zeros(1,8);
            res(1) = pFull(1);
            res(2) = pFull(4);
            res(3) = pFull(2);
            res(4) = pFull(5);
            res(5) = pFull(7);
            res(6) = pFull(12);
            res(7) = pFull(9);
            res(8) = pFull(14);
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
                  + obj.mix( diagParKE(obj.frag.Z(iatom),itype), ...
                  obj.fnar.KE(ran,ran), ...
                  obj.fdif.KE(ran,ran) );
               res(ran,ran) = res(ran,ran) - obj.frag.H1en(ran,ran,iatom) ...
                  + obj.mix( diagParEN(obj.frag.Z(iatom),itype), ...
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
         res = obj.KE(ienv);
         for iatom = 1:obj.natom
            res = res + obj.H1en(iatom);
         end
         if (ienv > 0)
            res = res + obj.frag.H1Env(:,:,ienv);
         end
      end
      function res = KE(obj,ienv)
         par = obj.par;
         diagParKE = zeros(6,2); % element and type(s,p)
         diagParKE(1,1) = par(1); % diagonal for H
         diagParKE(6,1) = par(2); % s diagonal for C
         diagParKE(6,2) = par(3); % p diagonal for C
         if (obj.rhodep == 1)
            diagParKErho = zeros(6,2);
            diagParKErho(1,1) = par(17);
            diagParKErho(6,1) = par(18);
            diagParKErho(6,2) = par(19);
         end
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
               mixpar = diagParKE(obj.frag.Z(iatom),itype);
               if (obj.rhodep == 1)
                   mixpar = mixpar + ...
                     diagParKErho(obj.frag.Z(iatom),itype) * ...
                     obj.rhozero(iatom,ienv+1);
               end
               ran = obj.valAtom{iatom,itype}; % range of valence orbs (1s)
               res(ran,ran) = res(ran,ran) - obj.frag.KE(ran,ran) ...
                  + obj.mix( mixpar, ...
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
               + obj.mix( diagParEN(obj.frag.Z(iatom),itype), ...
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
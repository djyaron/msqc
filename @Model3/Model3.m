classdef Model3 < handle
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
      aType   % (1,natom) atom type (initialized to Z)
      rcart   % (3,natom) cartesian coordinates of the atoms
      nenv
      X       % Transformation matrix used in hartreeFock
      
      nbasis  % number of atomic (and molecular) basis functions
      basisAtom  % (nbasis,1) atom # on which the function is centered
      basisType  % (nbasis,1) l quantum number: 0=s 1=p 2=d 3=d etc
      basisSubType % (nbasis,1) m quantum number: s=1 p=3 d=6 (cartesian)
      onAtom     % {natom,1}  list of basis functions on iatom
      valAtom    % {natom,type} list of valence basis functions of type
      %              (1-s 2-p) on iatom
      isBonded   % (natom,natom)  1 if atoms are bonded, 0 otherwise
      coord      % (natom,1)      # of atoms bonded to this atom
      charges    % (natom,nenv+1)  current charges on the atoms
      bondOrders % (natom,natom,nenv+1) current bond orders
      
      mixers     % (1,n)   cell array of mixers that are currently in use
      KEmods     % {1,n}   cell array of modifications to KE operator
      ENmods   % {natom}{1,n} cell array of modifications to H1en opers
      % a 1-elec operator modification has the following members
      %   ilist : modify ilist x jlist elements
      %   jlist :
      %   mixer  : pointer to a mix function
      
      % a 2-elec operator modification has the following members
      %   ilist, jlist, klist, llist:  elements to modify
      %   mixer : pointer to a mix function
      H2mods % {1,n}
      
      verify    % Verify output. Currently for C code.
      %end
      %properties (Transient)
      densitySave   % cell array {1:nenv+1} of most recent density matrices
      % used to start HF iterations
      % cached contexts (see atomContext and bondContext)
      atomContextXSaved % {iatom,ienv}
      atomContextNSaved % {iatom}
      bondContextXSaved % {iatom,jatom,ienv}
      bondContextNSaved % {iatom,jatom}
      index % used and managed externally
   end
   methods (Static)
      h2 = H2slater(F0, G1, F2)
   end
   methods
      function res = Model3(frag_,fnar_, fdif_)
         if (nargin ~= 0)
            res.frag = frag_;
            res.fnar = fnar_;
            res.fdif = fdif_;
            res.natom = frag_.natom;
            res.nelec = frag_.nelec;
            res.Z     = frag_.Z;
            res.aType = res.Z;
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
            res.coord = zeros(res.natom,1);
            for iatom = 1:res.natom
               for jatom = 1:res.natom
                  res.isBonded(iatom,jatom) = res.bonded(iatom,jatom);
               end
               res.coord(iatom) = sum(res.isBonded(iatom,:));
            end
            res.KEmods = cell(0,0);
            res.ENmods = cell(1,res.natom);
            for i=1:res.natom
               res.ENmods{1,i} = cell(0,0);
            end
            res.H2mods = cell(0,0);
            res.densitySave = cell(1,res.nenv+1);
            res.mixers = cell(0,0);
            if (isfield(frag_,'savedCharges'))
               res.charges = frag_.savedCharges;
            else
               % Initialize charges and bond orders
               res.charges = zeros(res.natom,res.nenv+1);
               for ienv = 0:res.nenv
                  res.charges(:,ienv+1) = res.frag.mcharge(ienv)';
               end
            end
            % Initialize bond orders
            if (isfield(frag_,'savedBondOrders'))
               res.bondOrders = frag_.savedBondOrders;
            else
               res.bondOrders = zeros(res.natom,res.natom,res.nenv+1);
               for ienv = 0:res.nenv
                  res.bondOrders(:,:,ienv+1) = res.frag.calcBO(ienv);
               end
            end
            res.EhfEnv  = zeros(1,res.nenv);
            res.EorbEnv = zeros(res.nbasis,res.nenv);
            res.orbEnv  = zeros(res.nbasis,res.nbasis,res.nenv);
            res.atomContextXSaved = {};
            res.atomContextNSaved = {};
            res.bondContextXSaved = {};
            res.bondContextNSaved = {};
         end
      end
      function clearModifiers(obj)
         obj.KEmods = cell(0,0);
         obj.ENmods = cell(1,obj.natom);
         for i=1:obj.natom
            obj.ENmods{1,i} = cell(0,0);
         end
         obj.H2mods = cell(0,0);
         obj.mixers = cell(0,0);
      end
      function res = npar(obj)
         res = 0;
         for i=1:size(obj.mixers,2)
            res = res + obj.mixers{1,i}.npar;
         end
      end
      function addMixer(obj, mix)
         add = 1;
         for i=1:size(obj.mixers,2)
            if (mix == obj.mixers{1,i})
               add = 0;
               break;
            end
         end
         if (add == 1)
            obj.mixers{1,end+1} = mix;
         end
      end
      function setPars(obj, pars)
         ic = 1;
         for i = 1:size(obj.mixers,2)
            mtemp = obj.mixers{1,i};
            n = mtemp.npar;
            if (n > 0)
               mtemp.setPars( pars(ic:(ic+n-1)));
            end
            ic = ic + n;
         end
      end
      function res = H1(obj, ienv)
         if (nargin < 2)
            ienv = 0;
         end
         res = obj.KE(ienv);
         for iatom = 1:obj.natom
            res = res + obj.H1en(iatom,ienv);
         end
         if (ienv > 0)
            res = res + obj.frag.H1Env(:,:,ienv);
         end
      end
      function res = KE(obj,ienv)
         % start with H1 matrix of unmodified STO-3G
         res   = obj.frag.KE;
         for imod = 1:size(obj.KEmods,2)
            mod = obj.KEmods{1,imod};
            ii = mod.ilist;
            jj = mod.jlist;
            tmp = mod.mixer.mix(obj.frag.KE(ii, jj), obj, ii, jj, ienv);
            res(ii,jj) = res(ii,jj) - obj.frag.KE(ii,jj) ...
               + tmp;
         end
      end
      %       function mixUsed = addKEmodConst(obj,mix)
      %           mod.ilist = 1:obj.nbasis;
      %           mod.jlist = 1:obj.nbasis;
      %           mod.mixer = mix;
      %           obj.KEmods{1,end+1} = mod;
      %           obj.addMixer(mix);
      %           mixUsed = mix;
      %       end
      function mixUsed = addKEmodDiag(obj,Zs,types,mix)
         if (nargin < 3)
            types = [1 2];
         end
         if (nargin < 4)
            mix = Mixer;
            % create a mix object for these blocks
            mix.desc = ['KE Diag Zs [',num2str(Zs),'] types [', ...
               num2str(types),']'];
         end
         mixerAdded = 0;
         for iZ = Zs % loop over all desired elements
            for iatom = find(obj.Z == iZ) % loop over atoms of this element
               ilist = []; % orbitals of "types" on this atom
               for itype = types
                  ilist = [ilist obj.valAtom{iatom,itype}'];
               end
               % Create a modifier for this block of the matrix
               mod.ilist = ilist;
               mod.jlist = ilist;
               mod.mixer = mix;
               obj.KEmods{1,end+1} = mod;
               mixerAdded = 1;
            end
         end
         if (mixerAdded)
            obj.addMixer(mix);
            mixUsed = mix;
         else
            mixUsed = [];
         end
      end
      function mixUsed = addKEcore(obj,Z,mix)
         mixerAdded = 0;
         for iatom = find(obj.Z == Z) % loop over atoms of this element
            if (length(obj.onAtom{iatom}) ~= 5)
               error('Model3:addKEcore not called on atom with 5 basis functions');
            end
            s1 = obj.onAtom{iatom}(1);
            mod.ilist = s1;
            mod.jlist = s1;
            mod.mixer = mix;
            obj.KEmods{1,end+1} = mod;
            mixerAdded = 1;
         end
         if (mixerAdded)
            obj.addMixer(mix);
            mixUsed = mix;
         else
            mixUsed = [];
         end
      end
      function mixUsed = addKEmodBonded(obj,Z1,Z2,types1,types2, mix)
         if (nargin < 4)
            types1 = [1 2];
         end
         if (nargin < 5)
            types2 = [1 2];
         end
         if (nargin < 6)
            mix = Mixer();
            mix.desc = ['KE bonded Z [',num2str(Z1),'] types [', ...
               num2str(types1),'] with Z [',num2str(Z2),'] types [', ...
               num2str(types2),']'];
         end
         mixerAdded = 0;
         for iatom = 1:obj.natom
            for jatom = 1:obj.natom
               bondExists = obj.isBonded(iatom,jatom);
               if ((mix.hybrid == 2) && bondExists)
                  % check for pi bond
                  bondExists = (obj.coord(iatom) == 3) && ...
                     (obj.coord(jatom) == 3);
               end
               if (bondExists)
                  for itype = 1:2
                     for jtype = 1:2
                        addmods = 0;
                        if ((obj.Z(iatom) == Z1) && (obj.Z(jatom) == Z2))
                           if (any(ismember(itype,types1)) && ...
                                 any(ismember(jtype,types2)) )
                              addmods = 1;
                           end
                        end
                        if ((obj.Z(iatom) == Z2) && (obj.Z(jatom) == Z1))
                           if (any(ismember(itype,types2)) && ...
                                 any(ismember(jtype,types1)) )
                              addmods = 1;
                           end
                        end
                        if (addmods)
                           mixerAdded = 1;
                           mod.ilist = obj.valAtom{iatom,itype}';
                           mod.jlist = obj.valAtom{jatom,jtype}';
                           mod.mixer = mix;
                           obj.KEmods{1,end+1} = mod;
                        end
                     end
                  end
               end
            end
         end
         if (mixerAdded)
            obj.addMixer(mix);
            mixUsed = mix;
         else
            mixUsed = [];
         end
      end
      function mixUsed = addMMStretch(obj,oper,Z1,Z2, mix)
         if (nargin < 5)
            mix = Mixer();
            mix.type = 2;
            mix.par = [0 0 0];
            mix.fixed = [0 0 0];
            mix.desc = ['MM ',oper,' Z [',num2str(Z1), ...
               '] with Z [',num2str(Z2),']'];
         end
         mixerAdded = 0;
         for iatom = 1:obj.natom
            for jatom = (iatom + 1):obj.natom
               bondExists = obj.isBonded(iatom,jatom);
               if (bondExists)
                  mixerAdded = 1;
                  mod.ilist = (1:obj.nbasis)';
                  mod.jlist = (1:obj.nbasis)';
                  mod.mixer = mix;
                  switch oper
                     case 'KE'
                        obj.KEmods{1,end+1} = mod;
                     otherwise
                        error('Model3.addMMStretch: unknown oper type');
                  end
               end
            end
         end
         if (mixerAdded)
            obj.addMixer(mix);
            mixUsed = mix;
         else
            mixUsed = [];
         end
      end
      function mixUsed = addKEmodBondedh(obj,Z1,Z2,mix)
         mixerAdded = 0;
         for iatom = 1:obj.natom
            for jatom = 1:obj.natom
               bondExists = obj.isBonded(iatom,jatom);
               if ((mix.hybrid == 2) && bondExists)
                  % check for pi bond
                  bondExists = (obj.coord(iatom) == 3) && ...
                     (obj.coord(jatom) == 3);
               end
               if (bondExists)
                  addmods = 0;
                  if ((obj.Z(iatom) == Z1) && (obj.Z(jatom) == Z2))
                     addmods = 1;
                  end
                  if ((obj.Z(iatom) == Z2) && (obj.Z(jatom) == Z1))
                     addmods = 1;
                  end
                  if (addmods)
                     mixerAdded = 1;
                     mod.ilist = [obj.valAtom{iatom,1}',...
                        obj.valAtom{iatom,2}'];
                     mod.jlist = [obj.valAtom{jatom,1}',...
                        obj.valAtom{jatom,2}'];
                     mod.mixer = mix;
                     obj.KEmods{1,end+1} = mod;
                  end
               end
            end
         end
         if (mixerAdded)
            obj.addMixer(mix);
            mixUsed = mix;
         else
            mixUsed = [];
         end
      end
      function res = H1en(obj, iatom, ienv)
         % start with H1 matrix of unmodified STO-3G
         res   = obj.frag.H1en(:,:,iatom);
         
         mods = obj.ENmods{1,iatom};
         for imod = 1:size(mods,2)
            mod = mods{1,imod};
            ii = mod.ilist;
            jj = mod.jlist;
            tmp = mod.mixer.mix(obj.frag.H1en(ii, jj, iatom), obj, ii, jj, ienv);
            res(ii,jj) = res(ii,jj) - obj.frag.H1en(ii,jj,iatom) ...
               + tmp;
         end
      end
      %       function mixUsed = addENmodConst(obj,mix)
      %          mixerAdded = 0;
      %          for iZ = Zs % loop over all desired elements
      %             for iatom = find(obj.Z == iZ) % loop over atoms of this element
      %                ilist = obj.onAtom{iatom}'; % orbitals on this atom
      %                % Create a modifier for this block of the matrix
      %                mod.ilist = ilist;
      %                mod.jlist = ilist;
      %                mod.mixer = mix;
      %                obj.ENmods{1,end+1} = mod;
      %                mixerAdded = 1;
      %             end
      %          end
      %          if (mixerAdded)
      %             mixUsed = mix;
      %             obj.addMixer(mix);
      %          else
      %             mixUsed = [];
      %          end
      %       end
      function mixUsed = addENmodDiag(obj,Zs,types,mix)
         if (nargin < 3)
            types = [1 2];
         end
         if (nargin < 4)
            mix = Mixer;
            mix.desc = ['EN Diag Zs [',num2str(Zs),'] types [', ...
               num2str(types),']'];
         end
         mixerAdded = 0;
         % create a mix object that will be the same for all these blocks
         for iZ = Zs % loop over all desired elements
            for iatom = find(obj.Z == iZ) % loop over atoms of this element
               ilist = []; % orbitals of "types" on this atom
               for itype = types
                  ilist = [ilist obj.valAtom{iatom,itype}'];
               end
               % Create a modifier for this block of the matrix
               mod.ilist = ilist;
               mod.jlist = ilist;
               mod.mixer = mix;
               obj.ENmods{iatom}{1,end+1} = mod;
               mixerAdded = 1;
            end
         end
         if (mixerAdded)
            obj.addMixer(mix);
            mixUsed = mix;
         else
            mixUsed = [];
         end
      end
      function mixUsed = addENcore(obj,Z,mix)
         mixerAdded = 0;
         for iatom = find(obj.Z == Z) % loop over atoms of this element
            if (length(obj.onAtom{iatom}) ~= 5)
               error('Model3:addENcore not called on atom with 5 basis functions');
            end
            s1 = obj.onAtom{iatom}(1);
            mod.ilist = s1;
            mod.jlist = s1;
            mod.mixer = mix;
            obj.ENmods{iatom}{1,end+1} = mod;
            mixerAdded = 1;
         end
         if (mixerAdded)
            obj.addMixer(mix);
            mixUsed = mix;
         else
            mixUsed = [];
         end
      end
      function mixUsed = addENmodBonded(obj,Z1,Z2,types1,types2, mix)
         if (nargin < 4)
            types1 = [1 2];
         end
         if (nargin < 5)
            types2 = [1 2];
         end
         if (nargin < 6)
            mix = Mixer();
            mix.desc = ['EN bonded Z ',num2str(Z1),' types [', ...
               num2str(types1),'] with Z ',num2str(Z2),' types [', ...
               num2str(types2),']'];
         end
         mixerAdded = 0;
         for operAtom = 1:obj.natom
            for iatom = 1:obj.natom
               for jatom = 1:obj.natom
                  bondExists = obj.isBonded(iatom,jatom);
                  if ((mix.hybrid == 2) && bondExists)
                     % check for pi bond
                     bondExists = (obj.coord(iatom) == 3) && ...
                        (obj.coord(jatom) == 3);
                  end
                  if (bondExists && ...
                        ((iatom==operAtom) || (jatom==operAtom)) )
                     for itype = 1:2
                        for jtype = 1:2
                           addmods = 0;
                           if ((obj.Z(iatom) == Z1) && (obj.Z(jatom) == Z2))
                              if (any(ismember(itype,types1)) && ...
                                    any(ismember(jtype,types2)) )
                                 addmods = 1;
                              end
                           end
                           if ((obj.Z(iatom) == Z2) && (obj.Z(jatom) == Z1))
                              if (any(ismember(itype,types2)) && ...
                                    any(ismember(jtype,types1)) )
                                 addmods = 1;
                              end
                           end
                           if (addmods)
                              mod.ilist = obj.valAtom{iatom,itype}';
                              mod.jlist = obj.valAtom{jatom,jtype}';
                              mod.mixer = mix;
                              obj.ENmods{1,operAtom}{1,end+1} = mod;
                              mixerAdded = 1;
                           end
                        end
                     end
                  end
               end
            end
         end
         if (mixerAdded)
            obj.addMixer(mix);
            mixUsed = mix;
         else
            mixUsed = [];
         end
      end
      function mixUsed = addENmodBonded1(obj,Z1,Z2,types1,types2, mix)
         % Modifies only the EN operator for atoms that match Z1
         if (nargin < 4)
            types1 = [1 2];
         end
         if (nargin < 5)
            types2 = [1 2];
         end
         if (nargin < 6)
            mix = Mixer();
            mix.desc = ['EN bonded(1 only) Z ',num2str(Z1),' types [', ...
               num2str(types1),'] with Z ',num2str(Z2),' types [', ...
               num2str(types2),']'];
         end
         mixerAdded = 0;
         for iatom = find(obj.Z == Z1)
            for jatom = find(obj.Z == Z2)
               bondExists = obj.isBonded(iatom,jatom);
               if ((mix.hybrid == 2) && bondExists)
                  % check for pi bond
                  bondExists = (obj.coord(iatom) == 3) && ...
                     (obj.coord(jatom) == 3);
               end
               if (bondExists)
                  for itype = 1:2
                     for jtype = 1:2
                        if (any(ismember(itype,types1)) && ...
                              any(ismember(jtype,types2)) )
                           mod.ilist = obj.valAtom{iatom,itype}';
                           mod.jlist = obj.valAtom{jatom,jtype}';
                           mod.mixer = mix;
                           obj.ENmods{1,iatom}{1,end+1} = mod;
                           mod.jlist = obj.valAtom{iatom,itype}';
                           mod.ilist = obj.valAtom{jatom,jtype}';
                           obj.ENmods{1,iatom}{1,end+1} = mod;
                           mod.mixer = mix;
                           mixerAdded = 1;
                        end
                     end
                  end
               end
            end
         end
         if (mixerAdded)
            obj.addMixer(mix);
            mixUsed = mix;
         else
            mixUsed = [];
         end
      end
      function mixUsed = addENmodBonded1h(obj,Z1,Z2,mix)
         mixerAdded = 0;
         for iatom = find(obj.Z == Z1)
            for jatom = find(obj.Z == Z2)
               bondExists = obj.isBonded(iatom,jatom);
               if ((mix.hybrid == 2) && bondExists)
                  % check for pi bond
                  bondExists = (obj.coord(iatom) == 3) && ...
                     (obj.coord(jatom) == 3);
               end
               if (bondExists)
                  iatomList = [obj.valAtom{iatom,1}', ...
                     obj.valAtom{iatom,2}'];
                  jatomList = [obj.valAtom{jatom,1}', ...
                     obj.valAtom{jatom,2}'];
                  mod.ilist = iatomList;
                  mod.jlist = jatomList;
                  mod.mixer = mix;
                  obj.ENmods{1,iatom}{1,end+1} = mod;
                  mod.jlist = iatomList;
                  mod.ilist = jatomList;
                  obj.ENmods{1,iatom}{1,end+1} = mod;
                  mod.mixer = mix;
                  mixerAdded = 1;
               end
            end
         end
         if (mixerAdded)
            obj.addMixer(mix);
            mixUsed = mix;
         else
            mixUsed = [];
         end
      end
      function res = Hnuc(obj,ienv)
         if (ienv == 0)
            res = obj.frag.Hnuc;
         else
            res = obj.frag.HnucEnv(ienv);
         end
      end
      function mixUsed = addH2modDiag(obj,Z,mix)
         if (nargin < 3)
            mix = Mixer;
            % create a mix object for these blocks
            mix.desc = ['H2 Diag Zs [',num2str(Zs),']'];
         end
         mixerAdded = 0;
         for iatom = find(obj.Z == Z) % loop over atoms of this element
            % ilist = obj.onAtom{iatom}'; % orbitals on this atom
            ilist = [obj.valAtom{iatom,1}',obj.valAtom{iatom,2}'];
            % Create a modifier for this block of the matrix
            mod.ilist = ilist;
            mod.jlist = ilist;
            mod.klist = ilist;
            mod.llist = ilist;
            mod.mixer = mix;
            obj.H2mods{1,end+1} = mod;
            mixerAdded = 1;
         end
         if (mixerAdded)
            obj.addMixer(mix);
            mixUsed = mix;
         else
            mixUsed = [];
         end
      end
      function mixUsed = addH2core(obj,Z,mix)
         mixerAdded = 0;
         for iatom = find(obj.Z == Z) % loop over atoms of this element
            if (length(obj.onAtom{iatom}) ~= 5)
               error('Model3:addKEcore not called on atom with 5 basis functions');
            end
            s1 = obj.onAtom{iatom}(1);
            mod.ilist = s1;
            mod.jlist = s1;
            mod.klist = s1;
            mod.llist = s1;
            mod.mixer = mix;
            obj.H2mods{1,end+1} = mod;
            mixerAdded = 1;
         end
         if (mixerAdded)
            obj.addMixer(mix);
            mixUsed = mix;
         else
            mixUsed = [];
         end
      end
      function addH2modSlater(obj,Z,mixF0,mixG1,mixF2)
         mixerAdded = 0;
         for iatom = find(obj.Z == Z) % loop over atoms of this element
            % ilist = obj.onAtom{iatom}'; % orbitals on this atom
            ilist = [obj.valAtom{iatom,1}',obj.valAtom{iatom,2}'];
            if (length(ilist) ~= 4)
               error('using H2slater for element without 4 basis funcs');
            end
            % Create a modifier for this block of the matrix
            mod.ilist = ilist;
            mod.F0mixer = mixF0;
            mod.G1mixer = mixG1;
            mod.F2mixer = mixF2;
            obj.H2mods{1,end+1} = mod;
            mixerAdded = 1;
         end
         if (mixerAdded)
            obj.addMixer(mixF0);
            obj.addMixer(mixG1);
            obj.addMixer(mixF2);
         end
      end
      function mixUsed = addH2modOffDiag(obj,Z1,Z2, mix)
         if (nargin < 4)
            mix = Mixer();
            mix.desc = ['KE bonded Z ',num2str(Z1),' with Z ', ...
               num2str(Z2)];
         end
         mixerAdded = 0;
         for iatom = 1:obj.natom
            for jatom = 1:obj.natom
               if (iatom ~= jatom)
                  bondExists = obj.isBonded(iatom,jatom);
                  if (bondExists == mix.bonded)
                     if ( ((obj.Z(iatom) == Z1) && (obj.Z(jatom) == Z2)) || ...
                           ((obj.Z(iatom) == Z2) && (obj.Z(jatom) == Z1)) )
                        mixerAdded = 1;
                        ilist = [obj.valAtom{iatom,1}',obj.valAtom{iatom,2}'];
                        jlist = [obj.valAtom{jatom,1}',obj.valAtom{jatom,2}'];
                        mod.ilist = ilist;
                        mod.jlist = ilist;
                        mod.klist = jlist;
                        mod.llist = jlist;
                        mod.mixer = mix;
                        obj.H2mods{1,end+1} = mod;
                     end
                  end
               end
            end
         end
         if (mixerAdded)
            obj.addMixer(mix);
            mixUsed = mix;
         else
            mixUsed = [];
         end
      end
      function res = H2(obj,ienv)
         if (nargin < 2)
            ienv = 0;
         end
         res = obj.frag.H2;
         for imod = 1:length(obj.H2mods)
            mod = obj.H2mods{imod};
            if (isfield(mod,'jlist'))
               i = mod.ilist;
               j = mod.jlist;
               k = mod.klist;
               l = mod.llist;
               tmp = mod.mixer.mix(obj.frag.H2(i, j, k, l), obj, i, k, ienv);
               res(i,j,k,l) = res(i,j,k,l) - obj.frag.H2(i,j,k,l) ...
                  + tmp;
            else
               % F0 = h2(s,s,s,s);  F2 = h2(px,py,px,py)*25/3;
               % G1 = h2(s,px,s,px)*3;
               i = mod.ilist;
               s = i(1); px = i(2); py = i(3);
               F0 = mod.F0mixer.mix(obj.frag.H2(s,s,s,s), ...
                  obj, i, i, ienv);
               G1 = mod.G1mixer.mix(obj.frag.H2(s,px,s,px), ...
                  obj, i, i, ienv)*3;
               F2 = mod.F2mixer.mix(obj.frag.H2(px,py,px,py), ...
                  obj, i, i, ienv)*25/3;
               res(i,i,i,i) = obj.H2slater(F0,G1,F2);
            end
         end
      end
      function res = S(obj)
         res = obj.frag.S;
      end
      function res = dataForParallel(obj,scaleOnly)
         % copy data needed for parallel HF to a non-handle object
         % if scaleOnly = true, it only copies the frag data
         % otherwise, copies frag fdif and fnar
         fr.natom = obj.natom;
         fr.nelec = obj.nelec;
         fr.Z     = obj.Z;
         fr.rcart = obj.rcart;
         fr.nenv = obj.nenv;
         fr.nbasis = obj.nbasis;
         fr.basisAtom = obj.basisAtom;
         fr.basisType = obj.basisType;
         fr.basisSubType = obj.basisSubType;
         fr.savedCharges  = obj.charges;
         fr.savedBondOrders = obj.bondOrders;
         fr.KE   = obj.frag.KE;
         fr.H1en = obj.frag.H1en;
         fr.H2   = obj.frag.H2;
         fr.H1Env = obj.frag.H1Env;
         res.frag = fr;
         fn.KE   = obj.fnar.KE;
         fn.H1en = obj.fnar.H1en;
         fn.H2   = obj.fnar.H2;
         fn.H1Env = obj.fnar.H1Env;
         res.fnar = fn;
         fd.KE   = obj.fdif.KE;
         fd.H1en = obj.fdif.H1en;
         fd.H2   = obj.fdif.H2;
         fd.H1Env = obj.fdif.H1Env;
         res.fdif = fd;
         mixes = cell(1,length(obj.mixers));
         for i = 1:length(obj.mixers)
            obj.mixers{i}.index = i; % for use below
            mixes{i} = obj.mixers{i}.constructionData;
         end
         res.mixers = mixes;
         kemods = cell(1,length(obj.KEmods));
         for i = 1:length(obj.KEmods)
            t1 = [];
            t1.ilist = obj.KEmods{i}.ilist;
            t1.jlist = obj.KEmods{i}.jlist;
            t1.mixNum = obj.KEmods{i}.mixer.index;
            kemods{i} = t1;
         end
         res.KEmods = kemods;
         enmods = cell(1,obj.natom);
         for iatom = 1:obj.natom
            t1 = cell(1,length(obj.ENmods{iatom}));
            for i = 1:length(obj.ENmods{iatom})
               t2 = [];
               t2.ilist = obj.ENmods{iatom}{i}.ilist;
               t2.jlist = obj.ENmods{iatom}{i}.jlist;
               t2.mixNum = obj.ENmods{iatom}{i}.mixer.index;
               t1{i} = t2;
            end
            enmods{iatom} = t1;
         end
         res.ENmods = enmods;
         h2mods = cell(1,length(obj.H2mods));
         for i = 1:length(obj.H2mods)
            t1 = [];
            t1.ilist = obj.H2mods{i}.ilist;
            t1.jlist = obj.H2mods{i}.jlist;
            t1.klist = obj.H2mods{i}.klist;
            t1.llist = obj.H2mods{i}.llist;
            t1.mixNum = obj.H2mods{i}.mixer.index;
            h2mods{i} = t1;
         end
         res.H2mods = h2mods;
      end
   end % methods
   methods (Static)
      function res = createFromData(dat)
         res = Model3(dat.frag,dat.fnar,dat.fdif);
         res.mixers = cell(1,length(dat.mixers));
         for i = 1:length(dat.mixers)
            res.mixers{i} = Mixer.createFromData(dat.mixers{i});
         end
         res.KEmods = cell(1,length(dat.KEmods));
         for i = 1:length(dat.KEmods)
            res.KEmods{i}.ilist = dat.KEmods{i}.ilist;
            res.KEmods{i}.jlist = dat.KEmods{i}.jlist;
            res.KEmods{i}.mixer = res.mixers{dat.KEmods{i}.mixNum};
         end
         res.ENmods = cell(1,res.natom);
         for iatom = 1:res.natom
            t1 = cell(1,length(dat.ENmods{iatom}));
            for i = 1:length(dat.ENmods{iatom})
               t1{i}.ilist = dat.ENmods{iatom}{i}.ilist;
               t1{i}.jlist = dat.ENmods{iatom}{i}.jlist;
               t1{i}.mixer = res.mixers{dat.ENmods{iatom}{i}.mixNum};
            end
            res.ENmods{iatom} = t1;
         end
         res.H2mods = cell(1,length(dat.H2mods));
         for i = 1:length(dat.H2mods)
            res.H2mods{i}.ilist = dat.H2mods{i}.ilist;
            res.H2mods{i}.jlist = dat.H2mods{i}.jlist;
            res.H2mods{i}.klist = dat.H2mods{i}.klist;
            res.H2mods{i}.llist = dat.H2mods{i}.llist;
            res.H2mods{i}.mixer = res.mixers{dat.H2mods{i}.mixNum};
         end
      end
   end
end %

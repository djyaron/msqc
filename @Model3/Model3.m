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
      rcart   % (3,natom) cartesian coordinates of the atoms
      nenv
      
      H2j     % {nbasis,nbasis} cell array of coulomb integrals
      H2k     % {nbasis,nbasis} cell array of exchange integrals
      
      nbasis  % number of atomic (and molecular) basis functions
      basisAtom  % (nbasis,1) atom # on which the function is centered
      basisType  % (nbasis,1) l quantum number: 0=s 1=p 2=d 3=d etc
      basisSubType % (nbasis,1) m quantum number: s=1 p=3 d=6 (cartesian)
      onAtom     % {natom,1}  list of basis functions on iatom
      valAtom    % {natom,type} list of valence basis functions of type
      %              (1-s 2-p) on iatom
      isBonded   % (natom,natom)  1 if atoms are bonded, 0 otherwise
      charges    % (natom,nenv+1)  current charges on the atoms
      
      mixers     % (1,n)   cell array of mixers that are currently in use
      KEmods     % {1,n}   cell array of modifications to KE operator
      ENmods   % {natom}{1,n} cell array of modifications to H1en opers
      % a modification has the following members
      %   ilist : modify ilist x jlist elements
      %   jlist :
      %   mixer  : pointer to a mix function
      
   end
   properties (Transient)
      densitySave   % cell array {1:nenv+1} of most recent density matrices
      % used to start HF iterations
   end
   methods
      function res = Model3(frag_,fnar_, fdif_)
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
         res.KEmods = cell(0,0);
         res.ENmods = cell(1,res.natom);
         for i=1:res.natom
            res.ENmods{1,i} = cell(0,0);
         end
         res.densitySave = cell(1,res.nenv+1);
         res.mixers = cell(0,0);
         res.charges = zeros(res.natom,res.nenv+1);
         for ienv = 0:res.nenv
            res.charges(:,ienv+1) = res.frag.mcharge(ienv)';
         end
         res.H2j = cell(res.nbasis,res.nbasis);
         res.H2k = cell(res.nbasis,res.nbasis);
         for i=1:res.nbasis
            for j=1:res.nbasis
               res.H2j{i,j} = squeeze(frag_.H2(i,j,:,:));
               res.H2k{i,j} = squeeze(frag_.H2(i,:,:,j));
            end
         end
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
      function setPar(obj, pars)
         ic = 1;
         for i = 1:size(obj.mixers,2)
            mtemp = obj.mixers{1,i};
            mtemp.setPars( pars(ic:(ic+mtemp.npar-1)));
            ic = ic + mtemp.npar;
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
         % start with H1 matrix of unmodified STO-3G
         res   = obj.frag.KE;
         
         for imod = 1:size(obj.KEmods,2)
            mod = obj.KEmods{1,imod};
            ii = mod.ilist;
            jj = mod.jlist;
            res(ii,jj) = res(ii,jj) - obj.frag.KE(ii,jj) ...
               + mod.mixer.mix(obj.fnar.KE(ii,jj), obj.fdif.KE(ii,jj), ...
               obj,ii,jj,ienv);
         end
      end
      function mixUsed = addKEmodDiag(obj,Zs,types,mix)
         if (nargin < 3)
            types = [1 2];
         end
         if (nargin < 4)
            mix = Mixer;
            mix.desc = ['KE Diag Zs [',num2str(Zs),'] types [', ...
               num2str(types),']'];
         end
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
               obj.KEmods{1,end+1} = mod;
            end
         end
         obj.addMixer(mix);
         mixUsed = mix;
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
         
         for iatom = 1:obj.natom
            for jatom = 1:obj.natom
               if (obj.isBonded(iatom,jatom))
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
                           obj.KEmods{1,end+1} = mod;
                        end
                     end
                  end
               end
            end
         end
         obj.addMixer(mix);
         mixUsed = mix;
      end
      function res = H1en(obj, iatom)
         % start with H1 matrix of unmodified STO-3G
         res   = obj.frag.H1en(:,:,iatom);
         
         mods = obj.ENmods{1,iatom};
         for imod = 1:size(mods,2)
            mod = mods{1,imod};
            ii = mod.ilist;
            jj = mod.jlist;
            res(ii,jj) = res(ii,jj) - obj.frag.H1en(ii,jj,iatom) ...
               + mod.mixer.mix(obj.fnar.H1en(ii,jj,iatom), ...
               obj.fdif.H1en(ii,jj,iatom) );
         end
      end
      function mixUsed = addENmodDiag(obj,Zs,types,mix)
         if (nargin < 3)
            types = [1 2];
         end
         if (nargin < 4)
            mix = Mixer;
            mix.desc = ['EN Diag Zs [',num2str(Zs),'] types [', ...
               num2str(types),']'];
         end
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
            end
         end
         obj.addMixer(mix);
         mixUsed = mix;
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
         for operAtom = 1:obj.natom
            for iatom = 1:obj.natom
               for jatom = 1:obj.natom
                  if (obj.isBonded(iatom,jatom) && ...
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
                           end
                        end
                     end
                  end
               end
            end
         end
         
%          for iatom = 1:obj.natom
%             for jatom = 1:obj.natom
%                 if (obj.isBonded(iatom,jatom))
%                   if ( ((obj.Z(iatom) == Z1) && (obj.Z(jatom) == Z2)) || ...
%                         ((obj.Z(iatom) == Z2) && (obj.Z(jatom) == Z1)))
%                      ilist = [];
%                      for itype = types1
%                         ilist = [ilist obj.valAtom{iatom,itype}'];
%                      end
%                      jlist = [];
%                      for jtype = types2
%                         jlist = [jlist obj.valAtom{jatom,jtype}'];
%                      end
%                      
%                      mod.ilist = ilist;
%                      mod.jlist = jlist;
%                      mod.mixer = mix;
%                      obj.ENmods{1,iatom}{1,end+1} = mod;
%                      mod.ilist = jlist;
%                      mod.jlist = ilist;
%                      mod.mixer = mix;
%                      obj.ENmods{1,iatom}{1,end+1} = mod;
%                   end
%                 end
%             end
%          end
         obj.addMixer(mix);
         mixUsed = mix;
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
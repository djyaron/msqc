classdef MBandit < handle
   properties
      % data used to construct the model
      frag;  % Fragment with regular STO-3G basis set
      fnar;  % Fragment with narrow  STO-3G basis set
      fdif;  % Fragment with diffuse STO-3G basis set
      ienv;  % environment being considered (0=no environment)
      
      % Bookkeeping information copied or drived from frag
      natom   % number of atoms in fragment
      nelec   % number of electrons in the fragment
      Z       % (1,natom) atomic numbers of the molecules
      rcart   % (3,natom) cartesian coordinates of the atoms
      
      nbasis  % number of atomic (and molecular) basis functions
      basisAtom  % (nbasis,1) atom # on which the function is centered
      basisType  % (nbasis,1) l quantum number: 0=s 1=p 2=d 3=d etc
      basisSubType % (nbasis,1) m quantum number: s=1 p=3 d=6 (cartesian)
      onAtom     % {natom,1}  list of basis functions on iatom
      valAtom    % {natom,type} list of valence basis functions of type
      
      nbond   % number of bonds
      bond   % (2,nbond) list of all bonds, with bond being between
      % (1,ibond) and (2,ibond) with (1,ibond)<(2,ibond)
      
      % Context variables
      charges  % charges on each atom, induced by environment
      bondOrders % bond orders retrieved from frag at object creation
      
      % storage of two-electron integrals in convenient format for
      % hartree fock algorithm
      H2j     % {nbasis,nbasis} cell array of coulomb integrals
      H2k     % {nbasis,nbasis} cell array of exchange integrals
      
      % one electron matrix elements
      % These are initialized to those in frag, and then modified
      KE  % (nbasis,nbasis) Modified KE operator
      H1en % (nbasis,nbasis,natom) Modified elec-nuc operators
      
      % Hartree fock output (only valid after call to solveHF). 
      orb
      Eorb
      Ehf
      
   end
   methods
      function res = MBandit(frag_,fnar_, fdif_,ienv)
         res.frag = frag_;
         res.fnar = fnar_;
         res.fdif = fdif_;
         % Bookkeeping info copied or derived from frag
         res.natom = frag_.natom;
         res.nelec = frag_.nelec;
         res.Z     = frag_.Z;
         res.rcart = frag_.rcart;
         res.nbasis = frag_.nbasis;
         res.basisAtom = frag_.basisAtom;
         res.basisType = frag_.basisType;
         res.basisSubType = frag_.basisSubType;
         res.ienv = ienv;
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
         ibond = 0;
         for iatom = 1:res.natom
            for jatom = (iatom+1):res.natom
               if res.bonded(iatom,jatom)
                  ibond = ibond + 1;
                  res.bond(1,ibond) = iatom;
                  res.bond(2,ibond) = jatom;
               end
            end
         end
         res.nbond = ibond;
         % initialize H2 data for HF algorithm
         res.H2j = cell(res.nbasis,res.nbasis);
         res.H2k = cell(res.nbasis,res.nbasis);
         for i=1:res.nbasis
            for j=1:res.nbasis
               res.H2j{i,j} = squeeze(frag_.H2(i,j,:,:));
               res.H2k{i,j} = squeeze(frag_.H2(i,:,:,j));
            end
         end
         % Initialize charges
         res.charges =  res.frag.mcharge(res.ienv) - res.frag.mcharge(0);
         % Initialize bond orders
         res.bondOrders = res.frag.calcBO(res.ienv);
         % initialize KE and H1en
         res.KE = res.frag.KE;
         res.H1en = res.frag.H1en;
      end
      function res = H1(obj)
         res = obj.frag.KE;
         for iatom = 1:res.natom
            res = res + obj.H1en(:,:,iatom);
         end
         if (obj.ienv > 0)
            res = res + obj.frag.H1Env(:,:,ienv);
         end
      end
      function res = mix(obj,v1,v2,x)
         res = ((1.0-x)/2.0) * v1 + ((1.0+x)/2.0) * v2;
      end
      function res = getBasis(obj,iatom,types)
         % Nx1 vector of basis functions of requested types on an atom
         % iatom = atom number
         % types = list of types of orbitals (1=s, 2=p)
         %         if types = [], then returns all orbitals
         if (isempty(types))
            res = obj.onAtom{iatom};
         else
            res = [];
            for itype = types
               res = [res ; obj.valAtom{iatom,itype}];
            end
         end
      end
      function atomModifierKE(obj,iatom,types,x)
         % iatom = atom to modify
         % types = [] includes all valence orbs on that atom
         %         otherwise, does list of included orbitals
         % mixer = mixer to use for the modification
         basis = obj.getBasis(iatom,types);
         v0 = obj.frag.KE(basis,basis);
         v1 = obj.fnar.KE(basis,basis);
         v2 = obj.fdif.KE(basis,basis);
         obj.KE(basis,basis) = obj.KE(basis,basis) - ...
            v0 + obj.mix(v1,v2,x);
      end
      function atomModifierEN(obj,iatom,types,x)
         % iatom = atom to modify
         % mixer = mixer to use for the modification
         basis = obj.getBasis(iatom,types);
         v0 = obj.frag.H1en(basis,basis,iatom);
         v1 = obj.fnar.H1en(basis,basis,iatom);
         v2 = obj.fdif.H1en(basis,basis,iatom);
         obj.H1en(basis,basis,iatom) = obj.H1en(basis,basis,iatom) - ...
            v0 + obj.mix(v1,v2,x);
      end
      function bondModifierKE(obj,iatom,itypes,jatom,jtypes,x)
         % iatom,jatom = atoms to modify
         % types = [] includes all valence orbs on that atom
         %         otherwise, does list of included orbitals
         % mixer = mixer to use for the modification
         ibasis = obj.getBasis(iatom,itypes);
         jbasis = obj.getBasis(jatom,jtypes);
         v0 = obj.frag.KE(ibasis,jbasis);
         v1 = obj.fnar.KE(ibasis,jbasis);
         v2 = obj.fdif.KE(ibasis,jbasis);
         obj.KE(ibasis,jbasis) = obj.KE(ibasis,jbasis) - ...
            v0 + obj.mix(v1,v2,x);
         v0 = obj.frag.KE(jbasis,ibasis);
         v1 = obj.fnar.KE(jbasis,ibasis);
         v2 = obj.fdif.KE(jbasis,ibasis);
         obj.KE(jbasis,ibasis) = obj.KE(jbasis,ibasis) - ...
            v0 + obj.mix(v1,v2,x);
      end
      function bondModifierEN(obj,iatom,itypes,jatom,jtypes,x)
         % iatom,itypes = EN operator to modify
         % jatom,jtypes = off diagonal terms to this atom
         % mixer = mixer to use for the modification
         ibasis = obj.getBasis(iatom,itypes);
         jbasis = obj.getBasis(jatom,jtypes);
         v0 = obj.frag.H1en(ibasis,jbasis,iatom);
         v1 = obj.fnar.H1en(ibasis,jbasis,iatom);
         v2 = obj.fdif.H1en(ibasis,jbasis,iatom);
         obj.H1en(ibasis,jbasis,iatom) = obj.H1en(ibasis,jbasis,iatom) - ...
            v0 + obj.mix(v1,v2,x);
         v0 = obj.frag.H1en(jbasis,ibasis,iatom);
         v1 = obj.fnar.H1en(jbasis,ibasis,iatom);
         v2 = obj.fdif.H1en(jbasis,ibasis,iatom);
         obj.H1en(jbasis,ibasis,iatom) = obj.H1en(jbasis,ibasis,iatom) - ...
            v0 + obj.mix(v1,v2,x);
      end
      function res = atomContext(obj,iatom)
         % Context of atom: currently Z and charge
         res.Z = obj.Z(iatom);
         res.charge = obj.charges(iatom);
      end
      function res = bondContext(obj,ibond)
         % Context of atom: current Zs and bond order
         atom1 = obj.bond(1,ibond);
         atom2 = obj.bond(2,ibond);
         res.Zs = [obj.Z(atom1),obj.Z(atom2)];
         res.bondorder = obj.bondOrders(atom1,atom2);
      end
      function status = solveHF(obj)
         [obj.orb,obj.Eorb,obj.Ehf,status] = obj.hartreeFock();
      end
   end
end
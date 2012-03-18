function [CorrE, MP2, Ehf, Eorb, orb, Nelectrons, Z, rcart, ...
         dipole, mulliken, ...
          atom, type, subtype, nprims, prims ] = readfchk(fid1)
% reads contents of a formatted checkpoint file from Gaussian
% Input
% fid1: file handle
% Output
% Ehf: total SCF energy of the molecule
% Eorb: (nbasis,1) molecular orbital energies
% orb: (nbasis,nbasis) molecular orbital coefficients
% Nelectrons: number of electrons in the molecule
% Z: (1,natom) atomic numbers of the molecules
% rcart: (3,natom) cartesian coordinates of the atoms
% dipole: (3,1) dipole moment of molecule
% mulliken: (1,natom) mulliken charge on the atoms
%
% The following are (nbasis,1) arrays with info on the atomic basis
% atom: atom # on which the function is centered
% type: l quantum number: 0=s 1=p 2=d 3=d etc
% subtype: 1..(2*type+1)
% nprims: number of primitives in this function
% prims: {nbasis,1} cell array of matrices of size (2,nprims)
% with (1,:) being contraction coefficients and
% (2,:) being primimitive exponents

%fid1 = fopen(['data',filesep,'temp.fch']);
t1 = textscan(fid1,'%s');
text = t1{1};
%fclose(fid1);

% Will have findText() issue errors for us as appropriate
issueErrors = true;

% Orbital energies
phrase = {'Alpha','Orbital','Energies'};
loc = findText(text,phrase, issueErrors);
% The fifth word after alpha is the number of energies
Nenergies = str2num( text{loc+5} );
Eorb = zeros(Nenergies,1);

for i=1:Nenergies
   Eorb(i) = str2double(text{loc + 5 + i});
end

% Number of electrons
phrase = {'Number','of','electrons'};
loc = findText(text,phrase);
%The fourth word after number is the number of electrons
Nelectrons = str2num(text{loc + 4});

% Orbital coefficients
phrase = {'Alpha','MO','coefficients'};
loc = findText(text,phrase);
% The fifth word after alpha is the number of energies
Nvalues = str2num(text{loc+5});
temp = zeros(Nvalues,1);

for i=1:Nvalues
   temp(i) = str2double(text{loc + 5 + i});
end

orb = reshape(temp,Nenergies,Nenergies);

% The hartree fock energy is after SCF Energy R
phrase = {'SCF','Energy','R'};
loc = findText(text,phrase);

Ehf = str2double(text{loc+3});

% MP2 Energy
phrase = {'MP2','Energy','R'};
loc = findText(text,phrase);

if loc == 0
    MP2 = Ehf;
else
MP2 = str2double(text{loc+3});
end

%Correlation Energy
 CorrE = MP2-Ehf;

% the atomic numbers (Z) are after 'Atomic numbers'
phrase = {'Atomic','numbers'};
loc = findText(text,phrase);
natom = str2num(text{loc+4});
Z = zeros(1,natom);
for iatom=1:natom
   Z(1,iatom) = str2num(text{loc+4+iatom});
end

% Cartesian coordinates
phrase = {'Current','cartesian','coordinates'};
loc = findText(text,phrase);
rcart = zeros(3,natom);
icurr = loc + 6;
for iatom=1:natom
   for ix= 1:3
   rcart(ix,iatom) = str2double(text{icurr});
   icurr = icurr+1;
   end
end
% Convert from Bohr radii to Angstroms
rcart = rcart / 1.889726124565062;

% Dipole moment
phrase = {'Dipole','Moment'};
loc = findText(text,phrase);
n1 = str2num(text{loc+4});
dipole = zeros(n1,1);
for i=1:n1
   dipole(i,1) = str2num(text{loc+4+i});
end

% Mulliken charges
phrase = {'Mulliken','Charges'};
loc = findText(text,phrase);
n1 = str2num(text{loc+4});
mulliken = zeros(1,n1);
for i=1:n1
   mulliken(1,i) = str2num(text{loc+4+i});
end

% Basis set information
% First, read in the data as it is defined in the fchk file
phrase = {'Shell','types'};
loc = findText(text,phrase);
nshells = str2num(text{loc+4});
shellTypes = zeros(nshells,1);
for i=1:nshells
   shellTypes(i,1) = str2num(text{loc+4+i});
end

phrase = {'Number','of','primitives','per','shell'};
loc = findText(text,phrase);
nshells = str2num(text{loc+7});
primsPerShell = zeros(nshells,1);
for i=1:nshells
   primsPerShell(i,1) = str2num(text{loc+7+i});
end

phrase = {'Shell','to','atom','map'};
loc = findText(text,phrase);
nshells = str2num(text{loc+6});
shellToAtom = zeros(nshells,1);
for i=1:nshells
   shellToAtom(i,1) = str2num(text{loc+6+i});
end

phrase = {'Primitive','exponents'};
loc = findText(text,phrase);
n1 = str2num(text{loc+4});
primExp = zeros(n1,1);
for i=1:n1
   primExp(i,1) = str2num(text{loc+4+i});
end

phrase = {'Contraction','coefficients'};
loc = findText(text,phrase);
if (size(loc,1) > 0)
   n1 = str2num(text{loc(1)+4});
   contCoef = zeros(n1,1);
   for i=1:n1
      contCoef(i,1) = str2num(text{loc(1)+4+i});
   end
end
if (size(loc,1) == 2)
   if (strcmp(text{loc(2)-1},'P(S=P)')~=1)
      error('readfchk.m: unsupported contraction in fchk file');
   end
   n1 = str2num(text{loc(2)+4});
   contCoef2 = zeros(n1,1);
   for i=1:n1
      contCoef2(i,1) = str2num(text{loc(2)+4+i});
   end
end
if ((size(loc,2) > 2) || (size(loc,2) < 1))
  error('readfchk.m: unsupported contraction in fchk file');
end

% Information on the basis set. See top of file for definition of these
% arrays
atom = zeros(Nenergies, 1);
type = zeros(Nenergies, 1);
subtype = zeros(Nenergies, 1);
nprims = zeros(Nenergies, 1);
% prims{ibasis) = (2,nprims) matrix
% (1,:) = cont coefficients; (2,:) = prim exponents
prims = cell(Nenergies,1);

nsubtypes = [1 3 6]; %  cartesian basis sets for s,p,d
ibasis = 0;
iprim = 0;
for ishell = 1:nshells
   % stype = 0(s) 1(p) 2(d) etc.
   % if stype < 0, then it means we have shells up to that stype here
   % i.e. stype = -1, means we have both an s and p shell with identical
   % contractions
   stypeFile = shellTypes(ishell);
   if (stypeFile >= 0)
       stype = stypeFile;
       primTemp = zeros(2, primsPerShell(ishell));
       for ip = 1:primsPerShell(ishell)
           iprim = iprim+1;
           primTemp(1,ip) = contCoef(iprim);
           primTemp(2,ip) = primExp(iprim);
       end
       for itemp = 1:nsubtypes(stype+1)
           ibasis = ibasis + 1;
           atom(ibasis) = shellToAtom(ishell);
           type(ibasis) = stype;
           subtype(ibasis) = itemp;
           nprims(ibasis) = primsPerShell(ishell);
           prims{ibasis} = primTemp;
       end
   else
       for stype = 0:abs(stypeFile)
           primTemp = zeros(2, primsPerShell(ishell));
           for ip = 1:primsPerShell(ishell)
               iprim = iprim+1;
               if (stype == 0)
                  primTemp(1,ip) = contCoef(iprim);
               else
                  primTemp(1,ip) = contCoef2(iprim);
               end
               primTemp(2,ip) = primExp(iprim);
           end
           % reset iprim, to keep count ok
           if (stype < abs(stypeFile))
               iprim = iprim - primsPerShell(ishell);
           end
           for itemp = 1:nsubtypes(stype+1)
               ibasis = ibasis + 1;
               atom(ibasis) = shellToAtom(ishell);
               type(ibasis) = stype;
               subtype(ibasis) = itemp;
               nprims(ibasis) = primsPerShell(ishell);
               prims{ibasis} = primTemp;
           end
       end
   end
end
if (ibasis ~= Nenergies)
   error('readchk.m: bookkeeping error in creation of atom()');
end
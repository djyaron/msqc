function [Eorb, orb, atom, Nelectrons, Ehf] = oldreadfchk(fid1)
% This is probably all nonoptimal, and I did it pretty much by trial and
% error. The following makes a cell array of strings, holding each word in
% the input file.
%fid1 = fopen('try8.fch');
t1 = textscan(fid1,'%s');
text = t1{1};


% the energies are after "Alpha Orbital Energies" 
% This returns the indices of all occurences of the word Alpha
% (stolen from: http://arstechnica.com/civis/viewtopic.php?f=20&t=296197)
alpha = find(ismember(text,'Alpha')==1);
[nfound,junk] = size(alpha);
success = 0;

for i=1:nfound
   if (strcmp(text(alpha(i)+1),'Orbital') && ...%if this, then text contains the energies
         strcmp(text(alpha(i)+2),'Energies'))
      success = alpha(i);
   end
end

if (success == 0) 
   error('readchk: Alpha Orbital Energies not found');
end

% The fifth word after alpha is the number of energies
Nenergies = str2num(cell2mat(text(success+5)));
Eorb = zeros(Nenergies,1);

for i=1:Nenergies
   Eorb(i) = str2double(cell2mat(text(success + 5 + i)));%Reads the energies into an array
end


%returns all occurrences of the word Shell
shell = find(ismember(text,'Shell')==1);
[nfound,junk] = size(shell);
shellTypes = 0;  % location of text "Shell types"
shellToAtom = 0; % location of text "Shell to atom map"
for i = 1:nfound
   if (strcmp(text(shell(i) + 1),'types'))
      shellTypes = shell(i);
   end
   if (strcmp(text(shell(i) + 1),'to')&&...
         strcmp(text(shell(i) + 2), 'atom') && ...
         strcmp(text(shell(i) + 3),'map'))
      shellToAtom = shell(i);
   end
end

if (shellTypes == 0)
    error ('readchk.m: Shell types not found');
end
if (shellToAtom == 0)
    error ('readchk.m: Shell to atom map not found');
end

nshells = str2num(cell2mat(text(shellTypes + 4)));
atom = zeros(Nenergies, 1);
ibasis = 0;
for ishell = 1:nshells
   % stype = 0(s) 1(p) 2(d) etc.
   stype = abs(str2num(cell2mat(text(shellTypes + 4 + ishell))));
   satom = abs(str2num(cell2mat(text(shellToAtom + 6+ ishell))));
   for itemp = 1:(2*stype + 1)
      ibasis = ibasis + 1;
      atom(ibasis) = satom;
   end
end
if (ibasis ~= Nenergies)
   error('readchk.m: bookkeeping error in creation of atom()');
end


%returns the indices of all occurrences of the word "Number"
number = find(ismember(text,'Number')==1);
[nfound, junk] = size(number);
success = 0;
for i=1:nfound
    if (strcmp(text(number(i) + 1), 'of') && strcmp(text(number(i)+2),'electrons'))
        success = number(i);
    end
end

if (success == 0)
    error ('readchk: Number of electrons not found');
end

%The fourth word after number is the number of electrons
Nelectrons = str2num(cell2mat(text(success + 4)));


%%
% the energies are after "Alpha Orbital Energies" 
% This returns the indices of all occurences of the word Alpha
% (stolen from: http://arstechnica.com/civis/viewtopic.php?f=20&t=296197)
alpha = find(ismember(text,'Alpha')==1);
[nfound,junk] = size(alpha);
success = 0;

for i=1:nfound
   if (strcmp(text(alpha(i)+1),'MO') && ...
         strcmp(text(alpha(i)+2),'coefficients'))
      success = alpha(i);
   end
end

if (success == 0)
   error('readchk: Alpha MO coefficients not found');
end

% The fifth word after alpha is the number of energies
Nvalues = str2num(cell2mat(text(success+5)));
temp = zeros(Nvalues,1);

for i=1:Nvalues
   temp(i) = str2double(cell2mat(text(success + 5 + i)));
end

orb = reshape(temp,Nenergies,Nenergies);

% the hartree fock energy is after SCF Energy R
scf = find(ismember(text,'SCF')==1);
nfound = size(scf,1);
success = 0;

for i=1:nfound
   if (strcmp(text(scf(i)+1),'Energy') && ...
         strcmp(text(scf(i)+2),'R'))
      success = scf(i);
   end
end

if (success == 0)
   error('readchk: SCF Energy R not found');
end
Ehf = str2double(cell2mat(text(success+3)));

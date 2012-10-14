function res = analyzeBond(mod,ienv,atom1,atom2)
% Determine variables that characterize a bond from iatom to jatom

charge1 = mod.charges(atom1, ienv+1);
charge2 = mod.charges(atom2, ienv+1);
bondOrder = mod.bondOrders(atom1,atom2,ienv+1);

res = [charge1, bondOrder, charge2];

if (0)
% rot = <non-hybrid|hybrid>
rot1 = getRotationSigma(mod,atom1,atom2);
rot2 = getRotationSigma(mod,atom2,atom1);

% rhoNM = <nonhybridN|rho|nonhybridM>
fullRho = mod.density(ienv);
basis1 = [mod.valAtom{atom1,1};mod.valAtom{atom1,2}];
basis2 = [mod.valAtom{atom2,1};mod.valAtom{atom2,2}];
rho11 = fullRho(basis1,basis1);
rho12 = fullRho(basis1,basis2);
rho22 = fullRho(basis2,basis2);
% rho = <hybridN|rho|hybridM>
%     = <hybridN|nonhybridN><nonhybridN|rho|nonhybridM><nonhybridM|hybridM>
res = [ rot1'*rho11*rot1, rot1'*rho12*rot2, rot2'*rho22*rot2];
end
end

function rot = getRotationSigma(mod,a1,a2)
% rot1 and rot2 is the rotation matrix to hybrid orbitals for the bond
% between atom1 and atom2. rot = <non-hybrid|hybrid>
%    the p orbital on atom a1 pointing along r12 (x') can be obtained by:
%       |x'> = sum_j |j><j|x'> = sum_j |j> r12_j
%     so the new p-orbital is just a linear combination of the old
%     orbitals, with the coefficients being the vector r12_j
%        <j|x'> = r12_j
%     including the s orbitals, we get.
%        rot = <non-hybrid|hybrid> = <j|a> = [c1, c2*r12]
%     for sp3 c1 = 1/2         c2 = sqrt(3)/2
%     for sp2 c1 = 1/sqrt(3)   c2 = sqrt(2/3)
%     for sp  c1 = 1/sqrt(2)   c2 = sqrt(2)
if (mod.Z(a1) == 1)
   % if hydrogen, no need to form hybrids, rot = 1x1 identity matrix
   rot = 1;
else
   coord = mod.coord(a1);
   r12 = mod.rcart(:,a2) - mod.rcart(:,a1);
   r12 = r12/norm(r12);
   if (coord == 4) % sp3 hybridization
      rot = [0.5; (sqrt(3)/2.0) * r12];
   elseif (coord == 3) % sp2 hybridization
      rot = [1.0/sqrt(3.0); sqrt(1.5) * r12];
   elseif (coord == 2)
      rot = (1.0/sqrt(2.0)) * [1; r12];
   else
      error('Mixer: hybrid with coord = 1 for a non-hydrogen atom');
   end
end
end

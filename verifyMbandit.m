clear classes;
dataDir = './datasets';
load([dataDir,'/ethaneDat.mat']);
tic
for ienv = 1:100

% For comparison to older code
m3 = Model3(LL{1,1},LL{2,1},LL{3,1});

kH = rand;
kC = rand;
kCH = rand;
kCC = rand;
m3.addKEmodDiag(1,1,Mixer(kH,1,'ke.H'));
m3.addKEmodDiag(6,[1 2],Mixer(kC,1,'ke.C'));
m3.addKEmodBonded(1,6,1,[1 2],Mixer(kCH,1,'keCH'));
m3.addKEmodBonded(6,6,[1 2],[1 2],Mixer(kCC,1,'keCC'));

eH = rand;
eC = rand;
eCH = rand;
eCC = rand;
m3.addENmodDiag(1,1,Mixer(eH,1,'en.H'));
m3.addENmodDiag(6,[1 2],Mixer(eC,1,'en.C'));
m3.addENmodBonded(1,6,1,[1 2],Mixer(eCH,1,'en.CH'));
m3.addENmodBonded(6,6,[1 2],[1 2],Mixer(eCC,1,'en.CC'));

%m3.solveHF(ienv);

%
mb = Mbandit(LL{1,1},LL{2,1},LL{3,1},ienv);
for iatom = 1:mb.natom
   con = mb.atomContext(iatom);
   if (con.Z == 1)
      mb.atomModifierKE(iatom,1,kH);
   elseif (con.Z == 6)
      mb.atomModifierKE(iatom,[1 2],kC);
   end
end

for ibond = 1:mb.nbond
   atom1 = mb.bond(1,ibond);
   atom2 = mb.bond(2,ibond);
   con = mb.bondContext(ibond);
   Z1 = con.Zs(1);
   Z2 = con.Zs(2);
   if ((Z1 == 1) && (Z2 == 6))
      mb.bondModifierKE(atom1,1,atom2,[1 2],kCH);
   elseif ((Z1 == 6) && (Z2 == 1))
      mb.bondModifierKE(atom1,[1 2],atom2,1,kCH);
   elseif ((Z1 == 6) && (Z2 == 6))
      mb.bondModifierKE(atom1,[1 2],atom2,[1 2],kCC);
   end
end

for iatom = 1:mb.natom
   con = mb.atomContext(iatom);
   if (con.Z == 1)
      mb.atomModifierKE(iatom,1,eH);
   elseif (con.Z == 6)
      mb.atomModifierKE(iatom,[1 2],eC);
   end
end

for ibond = 1:mb.nbond
   atom1 = mb.bond(1,ibond);
   atom2 = mb.bond(2,ibond);
   con = mb.bondContext(ibond);
   Z1 = con.Zs(1);
   Z2 = con.Zs(2);
   if ((Z1 == 1) && (Z2 == 6))
      mb.bondModifierEN(atom1,1    ,atom2,[1 2],eCH);
      mb.bondModifierEN(atom2,[1 2],atom1,1    ,eCH);
   elseif ((Z1 == 6) && (Z2 == 1))
      mb.bondModifierEN(atom1,[1 2],atom2,1    ,eCH);
      mb.bondModifierEN(atom2,1    ,atom1,[1 2],eCH);
   elseif ((Z1 == 6) && (Z2 == 6))
      mb.bondModifierEN(atom1,[1 2],atom2,[1 2],eCC);
      mb.bondModifierEN(atom2,[1 2],atom1,[1 2],eCC);
   end
end

diffKE = max(max(abs(m3.KE(ienv) - mb.KE)));
for iatom = 1:mb.natom
   diffEN(iatom) = max(max(abs(m3.H1en(iatom,ienv) - mb.H1en(:,:,iatom))));
end

end
toc
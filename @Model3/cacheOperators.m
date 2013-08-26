function cacheOperators(obj,ienv)

cmix = changedMixers(obj,ienv);
changed = sum(cmix);
%disp(['mixers changed ',num2str(cmix(:)')]);
if (changed == 0)
elseif (changed == 1)
   imix = find(cmix == 1);
   updateOperatorsOneMixer(obj,ienv,imix);
   setCachedPars(obj,ienv);
else
   updateOperatorsFromScratch(obj,ienv);
   setCachedPars(obj,ienv);
end
end

function updateOperatorsFromScratch(obj,ienv)
obj.cacheValid(ienv+1) = 0;
t1 = obj.KE(ienv);
obj.cachedKE{ienv+1} = t1;
obj.cachedH1{ienv+1} = t1;
for iatom = 1:obj.natom
   t1 = obj.H1en(iatom,ienv);
   obj.cachedEN{iatom,ienv+1} = t1;
   obj.cachedH1{ienv+1} = obj.cachedH1{ienv+1} + t1;
end
if (ienv > 0)
   obj.cachedH1{ienv+1} = obj.cachedH1{ienv+1} + obj.frag.H1Env(:,:,ienv);
end

obj.cachedH2{ienv+1} = obj.H2(ienv);
obj.cacheValid(ienv+1) = 1;
end

function updateOperatorsOneMixer(obj,ienv,imix)
obj.cacheValid(ienv+1) = 0;
t1 = KE1mix(obj,ienv,imix);
obj.cachedKE{ienv+1} = t1;
obj.cachedH1{ienv+1} = t1;
for iatom = 1:obj.natom
   t1 = H1en1mix(obj,iatom,ienv,imix);
   obj.cachedEN{iatom,ienv+1} = t1;
   obj.cachedH1{ienv+1} = obj.cachedH1{ienv+1} + t1;
end
if (ienv > 0)
   obj.cachedH1{ienv+1} = obj.cachedH1{ienv+1} + obj.frag.H1Env(:,:,ienv);
end
obj.cachedH2{ienv+1} = H21mix(obj,ienv,imix);
obj.cacheValid(ienv+1) = 1;
end


function res = changedMixers(obj,ienv)
% return a list of all mixers that have changed

nmix = length(obj.mixers);
if (( length(obj.cachedPars) < ienv+1) || isempty(obj.cachedPars{ienv+1}))
   res = ones(nmix,1);
else
   res = zeros(nmix,1);
   pcurr = obj.cachedPars{ienv+1};
   pnew = cell(nmix,1);
   for imix = 1:nmix
      pnew{imix} = obj.mixers{imix}.par;
      if (sum(abs(pnew{imix}-pcurr{imix})) > 1.0e-12)
         res(imix) = 1;
      end
   end
end
end

function setCachedPars(obj,ienv)
nmix = length(obj.mixers);
pnew = cell(nmix,1);
for imix = 1:nmix
   pnew{imix} = obj.mixers{imix}.par;
end
obj.cachedPars{ienv+1} = pnew;
end

function res = KE1mix(obj,ienv,imix)
res   = obj.cachedKE{ienv+1};
%mixIndex = obj.mixers{imix}.index;
for imod = 1:length(obj.KEmods)
   mod = obj.KEmods{imod};
   %if (mod.mixer.index == mixIndex)
   if (mod.mixer == obj.mixers{imix})
      ii = mod.ilist;
      jj = mod.jlist;
      newPars = mod.mixer.par;
      oldPars = obj.cachedPars{ienv+1}{imix};
      % remove effect of cached mixer
      mod.mixer.par = oldPars;
      tmp = mod.mixer.mix(obj.frag.KE(ii, jj), obj, ii, jj, ienv);
      oldDelta = tmp - obj.frag.KE(ii,jj);
      mod.mixer.par = newPars;
      tmp = mod.mixer.mix(obj.frag.KE(ii, jj), obj, ii, jj, ienv);
      newDelta = tmp - obj.frag.KE(ii,jj);
      res(ii,jj) = res(ii,jj) - oldDelta + newDelta;
   end
end
end

function res = H1en1mix(obj, iatom, ienv, imix)
res   = obj.cachedEN{iatom,ienv+1};
mixIndex = obj.mixers{imix}.index;
mods = obj.ENmods{1,iatom};
for imod = 1:size(mods,2)
   mod = mods{imod};
   if (mod.mixer == obj.mixers{imix})
      ii = mod.ilist;
      jj = mod.jlist;
      newPars = mod.mixer.par;
      oldPars = obj.cachedPars{ienv+1}{imix};
      % remove effect of cached mixer
      mod.mixer.par = oldPars;
      tmp = mod.mixer.mix(obj.frag.H1en(ii, jj, iatom), obj, ii, jj, ienv);
      oldDelta = tmp - obj.frag.H1en(ii, jj, iatom);
      mod.mixer.par = newPars;
      tmp = mod.mixer.mix(obj.frag.H1en(ii, jj, iatom), obj, ii, jj, ienv);
      newDelta = tmp - obj.frag.H1en(ii, jj, iatom);
      res(ii,jj) = res(ii,jj) - oldDelta + newDelta;
   end
end
end

function res = H21mix(obj,ienv,imix)
res = obj.cachedH2{ienv+1};
mixIndex = obj.mixers{imix}.index;
for imod = 1:length(obj.H2mods)
   mod = obj.H2mods{imod};
   if (isfield(mod,'jlist'))
      if (mod.mixer == obj.mixers{imix})
         i = mod.ilist;
         j = mod.jlist;
         k = mod.klist;
         l = mod.llist;
         newPars = mod.mixer.par;
         oldPars = obj.cachedPars{ienv+1}{imix};
         % remove effect of cached mixer
         mod.mixer.par = oldPars;
         tmp = mod.mixer.mix(obj.frag.H2(i, j, k, l), obj, i, k, ienv);
         oldDelta = tmp - obj.frag.H2(i, j, k, l);
         mod.mixer.par = newPars;
         tmp = mod.mixer.mix(obj.frag.H2(i, j, k, l), obj, i, k, ienv);
         newDelta = tmp - obj.frag.H2(i, j, k, l);         
         res(i,j,k,l) = res(i,j,k,l) - oldDelta + newDelta;
      end
   else
      % F0 = h2(s,s,s,s);  F2 = h2(px,py,px,py)*25/3;
      % G1 = h2(s,px,s,px)*3;
      if ((mod.F0mixer == obj.mixers{imix}) || ...
            (mod.G1mixer == obj.mixers{imix}) || ...
            (mod.F2mixer == obj.mixers{imix}) )
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
end

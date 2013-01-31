function makeMixInfo(obj,atypes)

obj.atomTypes = atypes;
obj.mixInfo = cell(0,0);

for ip = 1:length(obj.policy)
   pol = obj.policy{ip};
   % policy has fields oper, func, sp, iatom, jatom, context
   if (~isfield(pol,'jatom')) % diagonal modifier
      for atype = getAtomTypes(pol.iatom,atypes)
         switch pol.oper
            case 'KE'
               obj.mixInfo = [obj.mixInfo, makeKEDiag(pol,atype)];
            case 'EN'
               obj.mixInfo = [obj.mixInfo, makeENDiag(pol,atype)]; 
            case 'E2'
               obj.mixInfo = [obj.mixInfo, makeE2Diag(pol,atype)]; 
             case '*'
               obj.mixInfo = [obj.mixInfo, makeKEDiag(pol,atype)];
               obj.mixInfo = [obj.mixInfo, makeENDiag(pol,atype)]; 
               obj.mixInfo = [obj.mixInfo, makeE2Diag(pol,atype)]; 
         end
      end
   else % off diagonal modifier
      for atype = getAtomTypes(pol.iatom,atypes)
         for btype = getAtomTypes(pol.jatom,atypes)
            if (atype < btype)
               continue
            end
            switch pol.oper
               case {'KE','EN'}
                  obj.mixInfo = [obj.mixInfo, ...
                      makeH1Bond(pol,pol.oper,atype,btype)];
               case 'E2'
                  obj.mixInfo = [obj.mixInfo, makeE2Bond(pol,atype,btype)];
                case '*'
                  obj.mixInfo = [obj.mixInfo, ...
                      makeH1Bond(pol,'KE',atype,btype)];
                  obj.mixInfo = [obj.mixInfo, ...
                      makeH1Bond(pol,'EN',atype,btype)];
                  obj.mixInfo = [obj.mixInfo, makeE2Bond(pol,atype,btype)];
                  
            end
         end
      end
   end
end

% Create array of mixers
obj.mixer = cell(0,0);
for i = 1:length(obj.mixInfo)
   switch obj.mixInfo{i}.type
      case 'E2slater'
         obj.addMixer(obj.mixInfo{i}.mixerF0);
         obj.addMixer(obj.mixInfo{i}.mixerG1);
         obj.addMixer(obj.mixInfo{i}.mixerF2);
      otherwise
         obj.addMixer(obj.mixInfo{i}.mixer);
   end
end

end

function res = getAtomTypes(itype,atypes)
if (strcmp(itype,'*'))
   res = atypes;
else
   res = itype;
end
end

function res = makeMixer(pol,desc,hybrid)
if (nargin < 3)
   hybrid = 0;
end
res = MixerC(0, pol.func, hybrid);
res.fixed = 0;
res.desc = desc;
if (isfield(pol,'jatom'))
   res.isDiag = 0;
   res.bonded = ~pol.nonbond;
else
   res.isDiag = 1;
end
if (isfield(pol,'context'))
   res.context = pol.context;
   nc = length(parseContext(res.context));
   res.par = zeros(1,1+nc);
   res.fixed = [0 ones(1,nc)];
end
end

function res = parseContext(cstr)
res = cell(0,0);
remain = cstr;
while (~isempty(remain))
   [res{end+1}, remain] = strtok(remain, ' ');
end
end

function res = makeInfo(mixer,type, iatom, jatom)
res.type = type;
res.mixer = mixer;
res.iatom = iatom;
if (nargin == 4)
   res.jatom = jatom;
end
end

function res = makeKEDiag(pol,atype)
res = cell(0,0);
desc = ['KE atype ',num2str(atype)];
if (Context.atypeToZtype(atype) == 1)
   desc = [desc,' 1s'];
   res = makeInfo(makeMixer(pol,desc),'KEdiags',atype);
else
   switch pol.sp
      case 'separate'
         desc1 = [desc,' 2s']; 
         m1 = makeInfo(makeMixer(pol,desc1),'KEdiags',atype);
         desc2 = [desc,' 2p']; 
         m2 = makeInfo(makeMixer(pol,desc2),'KEdiagp',atype);
         res = {m1 m2};
      case 'combine'
         mix = makeMixer(pol,[desc,' 2sp']);
         m1 = makeInfo(mix,'KEdiags',atype);
         m2 = makeInfo(mix,'KEdiagp',atype);
         res = {m1 m2};
      case 'core'
         desc1 = [desc,' core']; 
         res = makeInfo(makeMixer(pol,desc1),'KEcore',atype);
      case 'sonly'
         desc1 = [desc,' 2s']; 
         res = makeInfo(makeMixer(pol,desc1),'KEdiags',atype);
      case 'ponly'
         desc2 = [desc,' 2p']; 
         res = makeInfo(makeMixer(pol,desc2),'KEdiagp',atype);
      case 'shift'
         mix = makeMixer(pol,['KE shift']);
         m1 = makeInfo(mix,'KEdiags',1);
         m2 = makeInfo(mix,'KEcore', 6);
         m3 = makeInfo(mix,'KEdiags',6);
         m4 = makeInfo(mix,'KEdiagp',6);
         res = {m1 m2 m3 m4};         
      otherwise
         error('sp policy not compatible with KEdiag');
   end
end
end

function res = makeENDiag(pol,atype)
res = cell(0,0);
desc = ['EN atypes ',num2str(atype)];
if (Context.atypeToZtype(atype) == 1)
   desc = [desc,' 1s'];
   res = makeInfo(makeMixer(pol,desc),'ENdiags',atype);
else
   switch pol.sp
      case 'separate'
         desc1 = [desc,' 2s']; 
         m1 = makeInfo(makeMixer(pol,desc1),'ENdiags',atype);
         desc2 = [desc,' 2p']; 
         m2 = makeInfo(makeMixer(pol,desc2),'ENdiagp',atype);
         res = {m1 m2};
      case 'combine'
         mix = makeMixer(pol,[desc,' 2sp']);
         m1 = makeInfo(mix,'ENdiags',atype);
         m2 = makeInfo(mix,'ENdiagp',atype);
         res = {m1 m2};
      case 'core'
         desc1 = [desc,' core']; 
         res = makeInfo(makeMixer(pol,desc1),'ENcore',atype);
      case 'sonly'
         desc1 = [desc,' 2s']; 
         res = makeInfo(makeMixer(pol,desc1),'ENdiags',atype);
      case 'ponly'
         desc2 = [desc,' 2p']; 
         res = makeInfo(makeMixer(pol,desc2),'ENdiagp',atype);
      case 'shift'
         mix = makeMixer(pol,'EN shift');
         m1 = makeInfo(mix,'ENcore', atype);
         m2 = makeInfo(mix,'ENdiags', atype);
         m3 = makeInfo(mix,'ENdiagp',atype);
         res = {m1 m2 m3};         
      otherwise
         error('sp policy not compatible with ENdiag');
   end
end
end

function res = makeE2Diag(pol,atype)
res = cell(0,0);
desc = ['E2 atype ',num2str(atype)];
switch pol.sp
   case 'core'
      desc = [desc,' core'];
      res = makeInfo(makeMixer(pol,desc),'E2core',atype);
   case 'slater'
      if (Context.atypeToZtype(atype) == 1)
         res = makeInfo(makeMixer(pol,desc),'E2diag',atype);
      else
         res = [];
         res.type = 'E2slater';
         res.iatom = atype;
         res.mixerF0 = makeMixer(pol,[desc,' F0']);
         res.mixerG1 = makeMixer(pol,[desc,' G1']);
         res.mixerF2 = makeMixer(pol,[desc,' F2']);
      end
   case 'shift'
      mix = makeMixer(pol,'E2 shift');
      m1 = makeInfo(mix,'E2diag', 1);
      m2 = makeInfo(mix,'E2diag', 6);
      res = {m1 m2};
   otherwise
      res = makeInfo(makeMixer(pol,desc),'E2diag',atype);
end
end

function res = makeH1Bond(pol,oper,atype,btype)
res = cell(0,0);
desc = [oper,' atypes ',num2str(atype),' ',num2str(btype)];
Za = Context.atypeToZtype(atype);
Zb = Context.atypeToZtype(btype);
if ((Za == 1) && (Za==1))
   desc = [desc,' ss'];
   res = makeInfo(makeMixer(pol,desc),[oper,'bondss'],atype,btype);
else
   switch pol.sp
      case 'combine'
         mix = makeMixer(pol,[desc,' 2sp']);
         res = {makeInfo(mix,[oper,'bondss'],atype,btype)};
         if (Za > 1)
            res{end+1} = makeInfo(mix,[oper,'bondps'],atype,btype);
         end
         if (Zb > 1)
            res{end+1} = makeInfo(mix,[oper,'bondsp'],atype,btype);
         end
         if ((Za > 1) && (Zb > 1))
            res{end+1} = makeInfo(mix,[oper,'bondpp'],atype,btype);
         end
      case 'separate'
         mixss = makeMixer(pol,[desc,'ss']);
         mixsp = makeMixer(pol,[desc,'sp']);
         mixpp = makeMixer(pol,[desc,'pp']);
         res = {makeInfo(mixss,[oper,'bondss'],atype,btype)};
         if (Za > 1)
            res{end+1} = makeInfo(mixsp,[oper,'bondps'],atype,btype);
         end
         if (Zb > 1)
            res{end+1} = makeInfo(mixsp,[oper,'bondsp'],atype,btype);
         end
         if ((Za > 1) && (Zb > 1))
            res{end+1} = makeInfo(mixpp,[oper,'bondpp'],atype,btype);
         end
       case 'hybrid'
         mixSigma = makeMixer(pol,[desc,' sig'],1);
         res = {makeInfo(mixSigma,[oper,'bondh'],atype,btype)};
         if (strcmpi(oper,'EN') && (atype ~= btype))
            res{end+1} = makeInfo(mixSigma,[oper,'bondh'],btype,atype);
         end
         if ((Za > 1) && (Zb > 1)) % could be a pi bond
            mixPi = makeMixer(pol,[desc,' pi'],2);
            res{end+1} = makeInfo(mixPi,[oper,'bondh'],atype,btype);
            if (strcmpi(oper,'EN') && (atype ~= btype))
               res{end+1} = makeInfo(mixPi,[oper,'bondh'],btype,atype);
            end
         end
      otherwise
         error(['sp policy not compatible with ','oper','Bond']);
   end
end

end

function res = makeE2Bond(pol,atype,btype)
res = cell(0,0);
desc = ['E2 atypes ',num2str(atype),' ',num2str(btype)];
res = makeInfo(makeMixer(pol,desc),'E2bond',atype,btype);
end

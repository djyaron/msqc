function res = printEDetails(obj, ofile)
if (nargin < 2)
   ofile = 1;
end
res = cell(0,0);

molName = {'CH4' 'C2H6' 'C2H4' 'C3H8' 'C4H6','H2'};

[err pnum etype] = obj.err(obj.getPars);

err =err*627.509;
eke = err(etype==1);
eH  = err(etype==11);
eC  = err(etype==16);
e2  = err(etype==2);
etot = err(etype==3);
eNoCost = err(etype>0);
tot = norm(eNoCost)/sqrt(length(eNoCost));
ke = norm(eke)/sqrt(length(eke));
H = norm(eH)/sqrt(length(eH));
C = norm(eC)/sqrt(length(eC));
r2 = norm(e2)/sqrt(length(e2));
rtot = norm(etot)/sqrt(length(etot));
x.ke = ke; x.H = H; x.C = C; x.e2=r2; x.etot = rtot; x.name = 'all';
res{1} = x;
fprintf(ofile,'avg %5.3f ke %5.3f H %5.3f C %5.3f E2 %5.3f tot %5.3f\n',...
   tot, ke, H, C, r2, rtot);
if (length(unique(pnum(etype>0))) > 1)
   for ip = unique(pnum)
      eavg = err(pnum==ip);
      eke = err(etype==1 & pnum==ip);
      eH  = err(etype==11 & pnum==ip);
      eC  = err(etype==16 & pnum==ip);
      e2  = err(etype==2 & pnum==ip);
      etot  = err(etype==3 & pnum==ip);
      avg = norm(eavg)/sqrt(length(eavg));
      ke = norm(eke)/sqrt(length(eke));
      H = norm(eH)/sqrt(length(eH));
      C = norm(eC)/sqrt(length(eC));
      r2 = norm(e2)/sqrt(length(e2));
      rtot = norm(etot)/sqrt(length(etot));
      if ( ((ip-790) > 0) && ((ip-790) < 7) )
         fprintf(ofile,'>> %s ',molName{ip-790});
         x.name = molName{ip-790};
      else
         fprintf(ofile,'>> %i ',ip);
         x.name = ['p',num2str(ip)];
      end
      x.ke = ke; x.H = H; x.C = C; x.e2=r2; x.etot = rtot; 
      res{end+1} = x;

      fprintf(ofile,'avg %5.3f ke %5.3f H %5.3f C %5.3f E2 %5.3f tot %5.3f\n',...
         tot, ke, H, C, r2, rtot);
   end
end

%fprintf(ofile,'\n');
end


